/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2023.
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
import static org.junit.Assert.*;

public class TestParseMolWithAtomLists {

    @Test
    public void withAtomList() throws Exception {
        String mol = "\n  ACCLDraw09152014282D\n\n" +
                " 10 10  3  0  0  0  0  0  0  0999 V2000\n" +
                "   12.5547   -9.5481    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   11.7570   -9.5481    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   12.9561  -10.2274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   11.3247   -8.7890    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   11.3170  -10.2119    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   12.5624  -10.9298    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   12.9561   -8.8044    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   13.7640  -10.2274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   13.7975   -8.8044    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   14.2040   -9.5223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  2  1  1  0  0  0  0\n" +
                "  3  1  2  0  0  0  0\n" +
                "  4  2  2  0  0  0  0\n" +
                "  5  2  1  0  0  0  0\n" +
                "  6  3  1  0  0  0  0\n" +
                "  7  1  1  0  0  0  0\n" +
                "  8  3  1  0  0  0  0\n" +
                "  9  7  2  0  0  0  0\n" +
                " 10  9  1  0  0  0  0\n" +
                " 10  8  2  0  0  0  0\n" +
                "  4 F    2   8   7\n" +
                "  5 F    2   7   8\n" +
                "  6 F    2   7   8\n" +
                "M  ALS   4  2 F O   N   \n" +
                "M  ALS   5  2 F N   O   \n" +
                "M  ALS   6  2 F N   O   \n" +
                "M  END";

        Chemical c = Chemical.parseMol(mol);
        assertEquals(10, c.getAtomCount());
        assertTrue(c.hasCoordinates());
    }
}
