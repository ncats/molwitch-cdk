package org.openscience.cdk.geometry.cip;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.cdk.CdkChemicalImpl;
import gov.nih.ncats.molwitch.cdk.writer.Mdl2000WriterFactory;
import gov.nih.ncats.molwitch.io.ChemFormat;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImpl;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.geometry.cip.CIPTool.CIP_CHIRALITY;
import org.openscience.cdk.geometry.cip.rules.CIPLigandRule;
import org.openscience.cdk.geometry.cip.rules.CIPLigandRule2;
import org.openscience.cdk.geometry.cip.rules.ISequenceSubRule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;

/**
 * This modified version of {@link CIPTool} is just present temporarily to correct
 * a CIP bug in CDK. CDK's default CIP rule works well, but has a bug in the recursive
 * part where, if there's a tie between 2 ligands based on atomic number/mass, it will
 * then get all sub-ligands and compare them. The issue is that each set of sub-ligands
 * is sorted by a mass+atomic number only priority, and each sub-ligand from the parent
 * ligand is compared only to its matched counterpart, not to all potentially higher 
 * priority ligands.
 *  
 *  <pre>
 * Consider this case below (4S)-2,3,4,5-tetramethyloctane:  
 *  
 *     3'    5'
 *     M     C
 *     |  4  | 
 * M-C-C-[C]-C-C-C-M
 *   | 3  :  5
 *   M    M
 *        4'
 * 
 * M = Methyl
 * C = Carbon
 * : = Dashed bond
 * # = Locant 
 * </pre>
 *                 
 * The 4 postion [C] is a stereo center with S configuration (as demonstrated
 * with the : dashed bond to the methyl group below). The carbon atoms at 3, 4
 * and 5 positions are tied for priority based on first pass CIP rules (atom
 * number and mass). However, 3 and 5 are quickly seen as higher priority than
 * 4' in the first tie-break as 4' has no sub-ligands. 3 and 5, however, both
 * have 2 sub-ligands. 3 has [3',2] and 5 has [5',6] as sub-ligands. The next 
 * step CDK CIP rules do is to ORDER these sub-ligands and then compare them.
 * In this case, we need to order [3',2] and [5',6] individually. For [3',2],
 * this is done by comparing only the ATOMIC NUMBER+MASS between 3' and 2,
 * and since they are both carbon atoms, this will result in a tie. This means
 * that if the starting order of [3',2] was [3',2] it will remain that way. If
 * the starting order is [2,3'] it will remain [2,3']. Since the input order
 * of the ligands is based on read-order (order of bonds in the bond table in
 * the CTAB, for example), this means the ordering of the sub-ligands is effectively
 * non-deterministic. 2 is certainly higher priority than 3' in full CIP rules,
 * and 6 is higher priorty than 5' as well. In the above cases there are 4 ways
 * the ligand ordering can go
 * <pre>
 * Possibility 1: both sub-ligands in right order
 *    3 vs 5 =>
 *    [2,3'] vs [6,5'] => 
 *    2 vs 6 => 2 is higher priority
 *    therefore 3 is higher priority (correct)
 * 
 * Possibility 2: 3 sub-ligands wrong order, 5 right order
 *    3 vs 5 =>
 *    [3',2] vs [6,5'] => 
 *    3' vs 6 => 6 is higher priority
 *    therefore 5 is higher priority (incorrect)
 * 
 * Possibility 3: 3 sub-ligands RIGHT order, 5 wrong order
 *    3 vs 5 =>
 *    [2,3'] vs [5',6] => 
 *    2 vs 5' => 2 is higher priority
 *    therefore 3 is higher priority (correct)
 * 
 * Possibility 4: both sub-ligands in WRONG order
 *    3 vs 5 =>
 *    [3',2] vs [5',6] => 
 *    3' vs 5' => TIE, go to next
 *    2 vs 6 => 2 is higher priority
 *    therefore 3 is higher priority (correct)
 * </pre>
 * <p>
 * Notice here that 3 of the 4 possibilities above will give the correct
 * CIP ordering. This is indeed typical. The bug only applies in circumstances
 * where 2 or more "principal" ligands are tied for CIP priority based on
 * atom number and atomic weight AND those 2 ligands both have sub-ligands 
 * (non-hydrogens) AND at least 2 sub-ligands of those ligands are tied for
 * atomic weight and atomic mass. Even still, the majority of the time the
 * criteria is met, CDK will still produce the correct CIP labeling.
 * </p>
 * <p>
 * This class corrects the bug by switching how the ordering is done to instead
 * call the full CIP ordering for sub-ligands. This is potentially computationally
 * more expensive. The actual code here is largely just copy/pasted methods needed
 * to perform the same methods as {@link CIPTool} and a reference to a modified {@link CIPLigandRule}
 * (called {@link CIPLigandRule2}).
 * 
 * </p>
 * 
 * 
 *   
 * @author tyler
 *
 */
public class CIPToolMod {

	public static boolean USE_NEW_CENTRES=true;

    public static void setUseNewCentres(boolean centresRule) {
        USE_NEW_CENTRES=centresRule;
    }
	
	private static ISequenceSubRule<ILigand> cipRule = new CIPLigandRule2();

    private final static int MAX_RINGS = 5;

    /**
	 * GSRS-MODIFIED: Temporary bug fix for {@link CIPTool#label(IAtomContainer)}
	 * 
     * @param container structure to label
     */
    public static void label(IAtomContainer container, CdkChemicalImpl chemical) {
    	//Experimental new labeller
        int fragmentCount = 0;
        Iterator<CdkChemicalImpl> iterator = chemical.connectedComponents();
        while(iterator.hasNext() ){
            fragmentCount++;
            CdkChemicalImpl fragment = iterator.next();
            System.out.printf("fragment %d has %d atoms and %d bonds\n", fragmentCount, fragment.getAtomCount(),
                    fragment.getBondCount());
        }
        if( fragmentCount == 0) {
            fragmentCount = 1;
        }

        int ringCount = getSizeOfLargestRingSystem(chemical);
        System.out.printf("got ringCount %d and fragmentCount %d\n", ringCount, fragmentCount);

        if(ringCount <= MAX_RINGS) {
    		com.simolecule.centres.CdkLabeller.label(container);
    		return;
    	}else {
            System.out.println("NOT using USE_NEW_CENTRES");
	        for (IStereoElement stereoElement : container.stereoElements()) {
	            if (stereoElement instanceof ITetrahedralChirality) {
	                ITetrahedralChirality tc = (ITetrahedralChirality) stereoElement;
	                tc.getChiralAtom().setProperty(CDKConstants.CIP_DESCRIPTOR, getCIPChirality(container, tc).toString());
	            } else if (stereoElement instanceof IDoubleBondStereochemistry) {
	                IDoubleBondStereochemistry dbs = (IDoubleBondStereochemistry) stereoElement;
	                dbs.getStereoBond()
	                        .setProperty(CDKConstants.CIP_DESCRIPTOR, getCIPChirality(container, dbs).toString());
	            }
	        }
    	}

    }
    

    /**
	 * GSRS-MODIFIED: Temporary bug fix for {@link CIPTool#getCIPChirality(IAtomContainer, ITetrahedralChirality)}
	 * 
     * @param container structure to label
     */
    public static CIP_CHIRALITY getCIPChirality(IAtomContainer container, ITetrahedralChirality stereoCenter) {

        // the LigancyFourChirality is kind of redundant but we keep for an
        // easy way to get the ILigands array
        LigancyFourChirality tmp = new LigancyFourChirality(container, stereoCenter);
        Stereo stereo = stereoCenter.getStereo();

        int parity = permParity(tmp.getLigands());

        if (parity == 0) return CIP_CHIRALITY.NONE;
        if (parity < 0) stereo = stereo.invert();

        if (stereo == Stereo.CLOCKWISE) return CIP_CHIRALITY.R;
        if (stereo == Stereo.ANTI_CLOCKWISE) return CIP_CHIRALITY.S;

        return CIP_CHIRALITY.NONE;
    }
    
    /**
  	 * GSRS-MODIFIED: Temporary bug fix for {@link CIPTool#getCIPChirality(IAtomContainer, IDoubleBondStereochemistry)}
  	 * 
     * @param container structure to label
     */    
    public static CIP_CHIRALITY getCIPChirality(IAtomContainer container, IDoubleBondStereochemistry stereoCenter) {

        IBond stereoBond = stereoCenter.getStereoBond();
        IBond leftBond = stereoCenter.getBonds()[0];
        IBond rightBond = stereoCenter.getBonds()[1];

        // the following variables are usd to label the atoms - makes things
        // a little more concise
        //
        // x       y       x
        //  \     /         \
        //   u = v    or     u = v
        //                        \
        //                         y
        //
        IAtom u = stereoBond.getBegin();
        IAtom v = stereoBond.getEnd();
        IAtom x = leftBond.getOther(u);
        IAtom y = rightBond.getOther(v);

        Conformation conformation = stereoCenter.getStereo();

        ILigand[] leftLigands = getLigands(u, container, v);
        ILigand[] rightLigands = getLigands(v, container, u);

        if (leftLigands.length > 2 || rightLigands.length > 2) return CIP_CHIRALITY.NONE;

        // invert if x/y aren't in the first position
        if (!leftLigands[0].getLigandAtom().equals(x)) conformation = conformation.invert();
        if (!rightLigands[0].getLigandAtom().equals(y)) conformation = conformation.invert();

        int p = permParity(leftLigands) * permParity(rightLigands);

        if (p == 0) return CIP_CHIRALITY.NONE;

        if (p < 0) conformation = conformation.invert();

        if (conformation == Conformation.TOGETHER) return CIP_CHIRALITY.Z;
        if (conformation == Conformation.OPPOSITE) return CIP_CHIRALITY.E;

        return CIP_CHIRALITY.NONE;
    }
   
  
    /**
     * GSRS-MODIFIED: Temporary bug fix 
     * 
     * @param atom
     * @param container
     * @param exclude
     * @return
     */
    private static ILigand[] getLigands(IAtom atom, IAtomContainer container, IAtom exclude) {

        List<IAtom> neighbors = container.getConnectedAtomsList(atom);

        ILigand[] ligands = new ILigand[neighbors.size() - 1];

        int i = 0;
        for (IAtom neighbor : neighbors) {
            if (!neighbor.equals(exclude)) ligands[i++] = new Ligand(container, new VisitedAtoms(), atom, neighbor);
        }

        return ligands;
    }

    
    
    

    
    /**
     * Obtain the permutation parity (-1,0,+1) to put the ligands in descending
     * order (highest first). A parity of 0 indicates two or more ligands were
     * equivalent.
     *
     * @param ligands the ligands to sort
     * @return parity, odd (-1), even (+1) or none (0)
     */
    private static int permParity(final ILigand[] ligands) {
    	
        // count the number of swaps made by insertion sort - if duplicates
        // are found the parity is 0
        int swaps = 0;

        for (int j = 1, hi = ligands.length; j < hi; j++) {
            ILigand ligand = ligands[j];
            int i = j - 1;
            int cmp = 0;
            while ((i >= 0) && (cmp = cipRule.compare(ligand, ligands[i])) > 0) {
                ligands[i + 1] = ligands[i--];
                swaps++;
            }
            if (cmp == 0) // identical entries
                return 0;
            ligands[i + 1] = ligand;
        }

        // odd (-1) or even (+1)
        return (swaps & 0x1) == 0x1 ? -1 : +1;
    
    }

    public static int getSizeOfLargestRingSystem(CdkChemicalImpl chemical) {
        CdkChemicalImpl copy = chemical.deepCopy();
        int maxRings = 0;
        for(int i = copy.getBondCount()-1; i >=0; i--) {
            Bond bond = copy.getBond(i);
            if(!bond.isInRing() ) {
                copy.removeBond(i);
            }
        }

        for(int i = copy.getAtomCount()-1; i >=0; i--) {
            Atom atom = copy.getAtom(i);
            if(atom.getBondCount() == 0) {
                System.out.printf("atom %d of symbol %s has no bonds and will be deleted\n", i, atom.getSymbol());
                copy.removeAtom(i);
            }
        }

        Iterator<CdkChemicalImpl> iterator = copy.connectedComponents();
        int fragmentCount=0;
        while(iterator.hasNext() ){
            CdkChemicalImpl fragment = iterator.next();
            int currentRingTotal= fragment.getBondCount() - fragment.getAtomCount() + 1;
            fragmentCount++;
            System.out.printf("fragment %d has %d atoms and %d bonds\n", fragmentCount, fragment.getAtomCount(),
                    fragment.getBondCount());
            if( currentRingTotal >maxRings) {
                maxRings = currentRingTotal;
            }
        }
        return maxRings;
    }

    private void writeMolTemp(CdkChemicalImpl chem) {
        Mdl2000WriterFactory writerFactory = new Mdl2000WriterFactory();
        Path temp = null;
        try {
            temp = Files.createTempFile("temp", ".mol");
            FileOutputStream stream = new FileOutputStream(temp.toFile());
            ChemFormat.MolFormatSpecification DEFAULT_MOL_SPEC = new ChemFormat.MolFormatSpecification();
            ChemicalWriterImpl writer= writerFactory.newInstance(stream, DEFAULT_MOL_SPEC);
            writer.write(chem);
            writer.close();
            System.out.printf("wrote file to %s\n", temp);

        } catch (IOException e) {
            System.err.println("Error writing molecule to file: " + e.getMessage());
        }
    }
}
