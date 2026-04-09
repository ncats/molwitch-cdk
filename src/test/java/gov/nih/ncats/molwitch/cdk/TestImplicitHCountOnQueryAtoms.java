package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Chemical;
import org.junit.Test;

import java.io.IOException;
import java.util.logging.Logger;

import static org.junit.Assert.*;
public class TestImplicitHCountOnQueryAtoms {
    private static Logger logger = Logger.getLogger("TestImplicitHCountOnQueryAtoms");
    @Test
    public void parseSmartsWithQueryAtoms() throws Exception{
        //use 'Si' in the well-connected atom because the CDK toolkit assigns H count
        //  more permissively
        Chemical c = Chemical.parse("[#7,#8]~C1=c2c3c(O[Si]([#6])(O)C3=O)cc(O)c2=C(O)\\C=C/1");
        for(Atom a : c.getAtoms()){
            logger.fine(String.format("atom %s, impl H count: %d ", a.getSymbol(), a.getImplicitHCount()));
            assertTrue(a.getImplicitHCount() >=0);
        }
    }

    @Test
    public void getMassFromAromaticMolFile() throws IOException {
        String mol ="\n" +
                "  CDK     09172016003D\n" +
                "\n" +
                " 21 24  0  0  0  0  0  0  0  0999 V2000\n" +
                "   -3.4541   -2.7023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -3.1200   -1.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -4.2200   -0.2200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -3.8800    1.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -4.9770    2.2631    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -2.4500    1.6800    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.3500    0.6600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    1.3100    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    2.8100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -0.7506    4.1096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.7502    4.1089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.3500    0.6600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.4500    1.6800    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.8800    1.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    4.2200   -0.2200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.1200   -1.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.6900   -0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.7500   -1.9700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.3984   -3.3226    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -0.7500   -1.9700    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.6900   -0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  1  0  0  0  0\n" +
                "  2  3  4  0  0  0  0\n" +
                "  3  4  4  0  0  0  0\n" +
                "  4  5  1  0  0  0  0\n" +
                "  4  6  4  0  0  0  0\n" +
                "  6  7  4  0  0  0  0\n" +
                "  7  8  1  0  0  0  0\n" +
                "  8  9  1  0  0  0  0\n" +
                "  9 10  1  0  0  0  0\n" +
                " 10 11  1  0  0  0  0\n" +
                "  9 11  1  0  0  0  0\n" +
                "  8 12  1  0  0  0  0\n" +
                " 12 13  4  0  0  0  0\n" +
                " 13 14  4  0  0  0  0\n" +
                " 14 15  4  0  0  0  0\n" +
                " 15 16  4  0  0  0  0\n" +
                " 16 17  4  0  0  0  0\n" +
                " 12 17  4  0  0  0  0\n" +
                " 17 18  1  0  0  0  0\n" +
                " 18 19  2  0  0  0  0\n" +
                " 18 20  1  0  0  0  0\n" +
                " 20 21  1  0  0  0  0\n" +
                "  2 21  4  0  0  0  0\n" +
                "  7 21  4  0  0  0  0\n" +
                "M  END";

        assertEquals(300, (int) Chemical.parseMol(mol).getMass());
    }
}
