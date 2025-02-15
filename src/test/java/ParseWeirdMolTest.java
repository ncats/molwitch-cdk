/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2025.
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
import org.junit.Test;

import java.io.IOException;
import static org.junit.Assert.*;

public class ParseWeirdMolTest {

    @Test
    public void testMolWithIdAsNumber() throws IOException {
        String mol="13879741\n" +
                "  Marvin  01132112512D          \n" +
                "\n" +
                " 10 10  0  0  0  0            999 V2000\n" +
                "   15.2630   -4.8951    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   16.6919   -5.7202    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   15.2630   -2.4201    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   13.8341   -3.2451    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   15.9775   -4.4826    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   15.9775   -3.6576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   16.6919   -4.8951    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   15.2630   -3.2451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   14.5485   -3.6576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   14.5485   -4.4826    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  5  1  0  0  0  0\n" +
                "  1 10  1  0  0  0  0\n" +
                "  2  7  1  0  0  0  0\n" +
                "  3  8  2  0  0  0  0\n" +
                "  4  9  1  0  0  0  0\n" +
                "  5  6  1  0  0  0  0\n" +
                "  5  7  1  0  0  0  0\n" +
                "  6  8  1  0  0  0  0\n" +
                "  8  9  1  0  0  0  0\n" +
                "  9 10  2  0  0  0  0\n" +
                "M  END";

        Chemical chem = Chemical.parse(mol);
        assertEquals(10, chem.getAtomCount());

    }

    @Test
    public void mol2() throws IOException{
        String mol = "599690\n" +
                "  Marvin  01132101202D          \n" +
                "\n" +
                " 14 15  0  0  0  0            999 V2000\n" +
                "    3.7911   -1.8061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.6412   -3.9284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.2947   -3.1849    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    4.4553   -3.8299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    4.6122   -3.0200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.8937   -2.6247    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.7467   -4.1639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.7408   -4.7473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    5.0784   -4.3707    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    5.3921   -2.7508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    5.8583   -4.1016    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    6.0152   -3.2917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    5.5490   -1.9409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    6.4813   -4.6424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  6  2  0  0  0  0\n" +
                "  2  3  1  0  0  0  0\n" +
                "  2  4  1  0  0  0  0\n" +
                "  2  7  1  0  0  0  0\n" +
                "  2  8  1  0  0  0  0\n" +
                "  3  6  1  0  0  0  0\n" +
                "  4  5  1  0  0  0  0\n" +
                "  4  9  2  0  0  0  0\n" +
                "  5  6  1  0  0  0  0\n" +
                "  5 10  2  0  0  0  0\n" +
                "  9 11  1  0  0  0  0\n" +
                " 10 12  1  0  0  0  0\n" +
                " 10 13  1  0  0  0  0\n" +
                " 11 12  2  0  0  0  0\n" +
                " 11 14  1  0  0  0  0\n" +
                "M  END";

        Chemical chem = Chemical.parse(mol);
        assertEquals(14, chem.getAtomCount());
    }
}
