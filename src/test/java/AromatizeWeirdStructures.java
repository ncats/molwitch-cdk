/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2020.
 *
 * This work is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License as published by the Free Software Foundation;
 * either version 2.1 of the License, or (at your option) any later version.
 *
 * This work is distributed in the hope that it will be useful, but without any warranty;
 * without even the implied warranty of merchantability or fitness for a particular purpose.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 *  if not, write to:
 *
 *  the Free Software Foundation, Inc.
 *  59 Temple Place, Suite 330
 *  Boston, MA 02111-1307 USA
 */

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.cdk.CdkBond;
import org.junit.Test;
import org.openscience.cdk.interfaces.IBond;

import java.util.List;
import java.util.stream.Collectors;
import static org.junit.Assert.*;

public class AromatizeWeirdStructures {

    @Test
    public void doesntErrorOut() throws Exception{
        Chemical c = Chemical.parse("O=C([C@H](Cc1c2ccccc2[n]c1)N)O");
        c.aromatize();
        List<IBond> unsetBonds = c.bonds().map(CdkBond::getIBondFor).filter(b-> b.getOrder() == IBond.Order.UNSET).collect(Collectors.toList());
        assertTrue(c.toSmiles().contains("[nH]"));
        assertTrue(unsetBonds.toString(), unsetBonds.isEmpty());
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
        System.out.println(key);
        assertNotNull(key);
    }
}
