package gov.nih.ncats.molwitch.cdk.search;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.cdk.CdkUtil;
import gov.nih.ncats.molwitch.search.MolSearcher;
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.Collections;
import java.util.Optional;

public class CdkMolSearcher implements MolSearcher {
    private final Pattern pattern;

    public CdkMolSearcher(String smartsQuery){
        try {
            CdkUtil.parseSmarts(smartsQuery);
            pattern = SmartsPattern.create(smartsQuery, CdkUtil.getChemObjectBuilder());
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
    public CdkMolSearcher(Chemical chemical) {
        IAtomContainer container = CdkUtil.toAtomContainer(chemical);
        IAtomContainer query = hasQueryAtomsOrBonds(container)
                ? CdkUtil.asQueryAtomContainer(container)
                : withoutStereo(container);

        normalizeAromaticBondAtoms(query);

        pattern = Pattern.findSubstructure(query);
    }

    private static IAtomContainer withoutStereo(IAtomContainer container) {
        try {
            IAtomContainer copy = container.clone();
            copy.setStereoElements(Collections.emptyList());
            return copy;
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException("Could not copy query molecule", e);
        }
    }

    private static boolean hasQueryAtomsOrBonds(IAtomContainer container){
        for(IAtom a : container.atoms()){

            IAtom aa= AtomRef.deref(a);
            if(aa instanceof IQueryAtom || aa instanceof IPseudoAtom){
                return true;
            }
        }

        for(IBond b : container.bonds()){
            IBond ib= BondRef.deref(b);
            if(ib instanceof IQueryBond){
                return true;
            }
        }
        return false;
    }

    @Override
    public Optional<int[]> search(Chemical targetChemical) {
        IAtomContainer target = CdkUtil.toAtomContainer(targetChemical);
        try {
            int[] mapping = pattern.match(target);
            if(mapping.length == 0){
                return Optional.empty();
            }
            return Optional.of(mapping);
        } catch (Throwable e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static void normalizeAromaticBondAtoms(IAtomContainer mol) {
        for (IBond bond : mol.bonds()) {
            if (bond.isAromatic()) {
                bond.getBegin().setIsAromatic(true);
                bond.getEnd().setIsAromatic(true);
            }
        }
    }
}
