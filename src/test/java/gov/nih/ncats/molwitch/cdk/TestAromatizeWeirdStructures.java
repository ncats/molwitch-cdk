package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Chemical;
import org.junit.Test;
import org.openscience.cdk.interfaces.IBond;

import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;


import static org.junit.Assert.*;

public class TestAromatizeWeirdStructures {
    private static Logger logger = Logger.getLogger("TestAromatizeWeirdStructures");

    @Test()
    public void doesntErrorOut() throws Exception{
        Chemical c = Chemical.parse("O=C([C@H](Cc1c2ccccc2[n]c1)N)O");
        c.aromatize();
        List<IBond> unsetBonds = c.bonds().map(CdkBond::getIBondFor).filter(b-> b.getOrder() == IBond.Order.UNSET).collect(Collectors.toList());
        logger.finest(String.format("c.toSmiles(): %s", c.toSmiles()));
        assertTrue(c.toSmiles().contains("[nH]"));
        assertTrue(unsetBonds.toString(), unsetBonds.isEmpty());
    }

    @Test
    public void doesNotErrorAsMolfile() throws Exception {
        String inputMol = "\n" +
                "  ACCLDraw02052620062D\n" +
                "\n" +
                " 15 16  0  0  1  0  0  0  0  0999 V2000\n" +
                "   14.5277   -8.9699    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   13.0536  -10.9113    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   11.4485  -10.3810    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   10.5323  -10.2847    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0\n" +
                "   10.3411   -9.3840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    9.5432   -8.9234    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    9.5432   -8.0022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   10.3411   -7.5415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   11.1389   -8.0022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   11.1389   -8.9234    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   11.8232   -9.5394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   12.7237   -9.3498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   13.3396  -10.0356    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
                "   14.2416   -9.8456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   14.8572  -10.5310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                " 14  1  1  0  0  0  0\n" +
                " 13  2  1  1  0  0  0\n" +
                "  3 11  2  0  0  0  0\n" +
                "  4  3  1  0  0  0  0\n" +
                "  5  4  1  0  0  0  0\n" +
                "  5 10  1  0  0  0  0\n" +
                "  6  5  2  0  0  0  0\n" +
                "  7  6  1  0  0  0  0\n" +
                "  8  7  2  0  0  0  0\n" +
                "  9  8  1  0  0  0  0\n" +
                " 10  9  2  0  0  0  0\n" +
                " 11 10  1  0  0  0  0\n" +
                " 12 11  1  0  0  0  0\n" +
                " 13 12  1  0  0  0  0\n" +
                " 14 13  1  0  0  0  0\n" +
                " 15 14  2  0  0  0  0\n" +
                "M  END\n";
        Chemical c = Chemical.parseMol(inputMol);
        c.aromatize();
        List<IBond> unsetBonds = c.bonds().map(CdkBond::getIBondFor).filter(b-> b.getOrder() == IBond.Order.UNSET).collect(Collectors.toList());
        assertTrue(c.toSmiles().contains("[nH]"));
        assertTrue(unsetBonds.isEmpty());
    }

    @Test
    public void testKekulize5MemberCopy() throws IOException {
        String m = "\n"
                + "   JSDraw208012111172D\n"
                + "\n"
                + "  7  7  0  0  0  0              0 V2000\n"
                + "   18.5975   -7.8074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                + "   19.7965   -6.8094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                + "   21.2602   -7.3486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                + "   21.6844   -8.8500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                + "   23.2432   -8.9105    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
                + "   22.5570   -6.4815    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
                + "   23.7825   -7.4466    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                + "  1  2  1  0  0  0  0\n"
                + "  2  3  1  0  0  0  0\n"
                + "  3  4  2  0  0  0  0\n"
                + "  4  5  1  0  0  0  0\n"
                + "  3  6  1  0  0  0  0\n"
                + "  6  7  2  0  0  0  0\n"
                + "  5  7  1  0  0  0  0\n"
                + "M  END\n"
                + "";
        Chemical chem = Chemical.parse(m);
        chem.aromatize();
        chem=chem.copy();
        String molb = chem.toMol();
        assertTrue("Aromatized 5 membered ring should have aromatic bonds",
                molb.contains("  3  4  4  0  0  0  0\n"
                        + "  4  5  4  0  0  0  0\n"
                        + "  3  6  4  0  0  0  0\n"
                        + "  6  7  4  0  0  0  0\n"
                        + "  5  7  4  0  0  0  0")
                );
        chem.kekulize();
        String mol = chem.toMol();
        assertTrue("Kekulized 5 membered ring should have single/double bonds",
                mol.contains("  1  2  1  0  0  0  0\n"
                        + "  2  3  1  0  0  0  0\n"
                        + "  3  4  2  0  0  0  0\n"
                        + "  4  5  1  0  0  0  0\n"
                        + "  3  6  1  0  0  0  0\n"
                        + "  6  7  2  0  0  0  0\n"
                        + "  5  7  1  0  0  0  0")
                );
    }
    
    @Test
    public void parseMolWithUNIIAsName() throws Exception{
        String mol = "7PDD38B23D\n" +
                "  Marvin  08311509412D          \n" +
                "\n" +
                "  5  4  0  0  0  0            999 V2000\n" +
                "   14.8634   -8.7164    0.0000 He  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   14.1436   -9.5357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   13.4106   -9.3452    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   12.7552   -9.7129    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   12.0275   -9.5274    0.0000 He  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  2  1  1  0  0  0  0\n" +
                "  3  2  1  0  0  0  0\n" +
                "  4  3  1  0  0  0  0\n" +
                "  5  4  1  0  0  0  0\n" +
                "A    1\n" +
                "_R5\n" +
                "A    5\n" +
                "_R6\n" +
                "M  ISO  2   1   5   5   6\n" +
                "M  SPL  1   2   1\n" +
                "M  END";

        Chemical chem = Chemical.parse(mol);
        assertEquals(5, chem.getAtomCount());
    }

    @Test
    public void parseSymyxMol() throws Exception{
        String mol = "\n" +
                "  Symyx   08281518552D 1   1.00000     0.00000     0\n" +
                "\n" +
                " 16 13  0     0  0            999 V2000\n" +
                "   15.2500   -8.0313    0.0000 Al  0  1  0  0  0  0           0  0  0\n" +
                "   11.7350  -10.3495    0.0000 O   0  0  0  0  0  0           0  0  0\n" +
                "   11.7350   -9.1684    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "   10.7122   -8.5779    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "   10.7122   -7.3968    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "    9.6893   -6.8062    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "    8.6665   -7.3968    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "    8.6665   -8.5779    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "    9.6893   -9.1684    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "   11.7350   -6.8062    0.0000 O   0  0  0  0  0  0           0  0  0\n" +
                "   11.7350   -5.6251    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "   12.7579   -5.0346    0.0000 C   0  0  0  0  0  0           0  0  0\n" +
                "   10.7122   -5.0346    0.0000 O   0  0  0  0  0  0           0  0  0\n" +
                "   12.7579   -8.5779    0.0000 O   0  5  0  0  0  0           0  0  0\n" +
                "   17.5625   -7.9688    0.0000 O   0  5  0  0  0  0           0  0  0\n" +
                "   17.5625   -7.9688    0.0000 O   0  5  0  0  0  0           0  0  0\n" +
                "  2  3  2  0     0  0\n" +
                "  3  4  1  0     0  0\n" +
                "  4  5  2  0     0  0\n" +
                "  5  6  1  0     0  0\n" +
                "  6  7  2  0     0  0\n" +
                "  8  7  1  0     0  0\n" +
                "  4  9  1  0     0  0\n" +
                "  9  8  2  0     0  0\n" +
                "  5 10  1  0     0  0\n" +
                " 10 11  1  0     0  0\n" +
                " 11 12  1  0     0  0\n" +
                " 13 11  2  0     0  0\n" +
                "  3 14  1  0     0  0\n" +
                "M  CHG  4   1   3  14  -1  15  -1  16  -1\n" +
                "M  END";
        Chemical chem = Chemical.parse(mol);
        assertEquals(16, chem.getAtomCount());
    }

    @Test
    public void computeInchiFromMolWithRGroups() throws Exception{
        String mol= "\n\n" +
                "\n" +
                "  5  4  0  0  0  0  0  0  0  0999 V2000\n" +
                "    3.8399    0.2547    0.0000 _R2-2  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.4092    0.2505    0.0000 _R1-3  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.3362    0.6681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.0501    0.2547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.7641    0.6722    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  5  1  1  0  0  0  0\n" +
                "  3  4  1  0  0  0  0\n" +
                "  2  3  1  0  0  0  0\n" +
                "  4  5  1  0  0  0  0\n" +
                "M  END";

        String key = Chemical.parseMol(mol).toInchi().getKey();
        logger.fine(key);
        assertNotNull(key);
    }
}
