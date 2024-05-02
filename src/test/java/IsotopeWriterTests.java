/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2024.
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
import gov.nih.ncats.molwitch.inchi.InChiResult;
import gov.nih.ncats.molwitch.inchi.Inchi;
import org.junit.Test;

import static org.junit.Assert.*;

public class IsotopeWriterTests {

    @Test
    public void writeSmilesWithIsotope() throws Exception{
        String mm="\n" +
                "  Marvin  01132101482D          \n\n" +
                "  1  0  0  0  0  0            999 V2000\n" +
                "   11.1250   -3.7188    0.0000 I   0  5  0  0  0  0  0  0  0  0  0  0\n" +
                "M  CHG  1   1  -1\n" +
                "M  ISO  1   1 121\n" +
                "M  END";
        Chemical cc= Chemical.parse(mm);
        //prints out [I-] instead of [121I-]
        String actualSmiles = cc.toSmiles();
        assertTrue(actualSmiles, actualSmiles.contains("121I-"));
    }

    @Test
    public void writeInchiWithIsotope() throws Exception{
        String mm="\n" +
                "  Marvin  01132101482D          \n\n" +
                "  1  0  0  0  0  0            999 V2000\n" +
                "   11.1250   -3.7188    0.0000 I   0  5  0  0  0  0  0  0  0  0  0  0\n" +
                "M  CHG  1   1  -1\n" +
                "M  ISO  1   1 121\n" +
                "M  END";
        Chemical cc=Chemical.parse(mm);
        //correctly has 121 ISO
        System.out.println(cc.toMol());
        InChiResult inChiResult = cc.toInchi();
        String inchi = inChiResult.getInchi();
//        System.out.println(inchi);
        Chemical icopy= Inchi.toChemical(inchi);
        //+"\n"+ inChiResult.getAuxInfo()
        //does not have 121 ISO, has 0 ISO
        assertTrue("Molfile should have 121 ISO after inchi roundtrip", icopy.toMol().contains("ISO  1   1 121"));
//        System.out.println(icopy.toMol());
    }
}
