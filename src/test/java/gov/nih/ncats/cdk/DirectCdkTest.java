package gov.nih.ncats.cdk;

import org.freehep.graphicsbase.util.Assert;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.geometry.cip.CIPTool;
import org.openscience.cdk.geometry.cip.CIPToolMod;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.Map;

public class DirectCdkTest {

    @Test
    public void testStereoPerception() {
        String testSmiles = "CC[C@@H](C)[C@H](C)[C@@H](C)CC";
        int expectedAtomTotal = 10;
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser       smipar  = new SmilesParser(builder);
        SmilesGenerator    smigen  = new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.Stereo);
        try {
            IAtomContainer mol = smipar.parseSmiles(testSmiles);
            CIPToolMod.label(mol);
            for(int a =0; a < mol.getAtomCount(); a++) {
                System.out.printf("atom %d = %s stereoparity: %s\n", a, mol.getAtom(1), mol.getAtom(a).getStereoParity());
                System.out.println("properties:");
                if(mol.getAtom(a).getProperties().size() == 0) {
                    System.out.println("no properties!");
                }
                for(Map.Entry<Object, Object> entry : mol.getAtom(a).getProperties().entrySet() ) {
                    System.out.printf("key: %s = %s\n", entry.getKey(), entry.getValue());
                }
            }
            System.out.println(smigen.create(mol));
            Assert.assertEquals(expectedAtomTotal, mol.getAtomCount());
        } catch (InvalidSmilesException ex) {
            System.err.println("BAD SMILES: " + testSmiles);
        } catch (CDKException e) {
            System.err.println("CDK Exception: " + e.getMessage());
        }
    }
}
