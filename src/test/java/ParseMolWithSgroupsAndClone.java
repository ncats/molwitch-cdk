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
import gov.nih.ncats.molwitch.cdk.CdkUtil;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.sgroup.SgroupKey;

import java.util.Collection;
import static org.junit.Assert.*;
public class ParseMolWithSgroupsAndClone {

    @Test
    public void doesntThrowCastException() throws Exception{
        String mol= "\n\n  JSDraw207062014492D\n" +
                " 19 19  0  0  0  0              0 V2000\n" +
                "   34.7524  -11.8314    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   34.7524  -10.3407    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   33.4688   -9.6284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   32.2293  -10.3407    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   30.8822   -9.5130    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   29.6100  -10.2950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   28.2993   -9.2792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.6430   -9.2792    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   27.5017   -7.8487    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   25.4978   -8.5482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   25.8524  -10.4985    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
                "   20.7800  -17.9621    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   32.2293  -11.8314    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   33.5320  -12.5242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   33.4688   -5.5321    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   36.0925   -9.4672    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   36.0925   -7.9686    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   36.0925  -12.6299    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   41.4798  -17.1579    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  2  0  0  0  0\n" +
                "  2  3  1  0  0  0  0\n" +
                "  3  4  2  0  0  0  0\n" +
                "  4  5  1  0  0  0  0\n" +
                "  5  6  1  0  0  0  0\n" +
                "  6  7  1  0  0  0  0\n" +
                "  7  8  1  0  0  0  0\n" +
                "  8  9  2  0  0  0  0\n" +
                "  8 10  2  0  0  0  0\n" +
                "  8 11  1  0  0  0  0\n" +
                "  6 12  1  0  0  0  0\n" +
                " 13  4  1  0  0  0  0\n" +
                " 14 13  2  0  0  0  0\n" +
                "  1 14  1  0  0  0  0\n" +
                "  3 15  1  0  0  0  0\n" +
                "  2 16  1  0  0  0  0\n" +
                " 16 17  1  0  0  0  0\n" +
                "  1 18  1  0  0  0  0\n" +
                " 18 19  1  0  0  0  0\n" +
                "M  STY  1   1 SRU\n" +
                "M  SCN  1   1 HT \n" +
                "M  SAL   1  8   1   2   3   4   5   6   7   8\n" +
                "M  SAL   1  8   9  10  11  13  14  15  16  17\n" +
                "M  SAL   1  1  18\n" +
                "M  SBL   1  2  11  19\n" +
                "M  SMT   1 A\n" +
                "M  SPA   1  8   1   2   3   4   5   6   7   8\n" +
                "M  SPA   1  8   9  10  11  13  14  15  16  17\n" +
                "M  SPA   1  1  18\n" +
                "M  SDI   1  4   25.1727  -16.0952   25.1727  -12.8501\n" +
                "M  SDI   1  4   38.7614  -12.8501   38.7614  -16.0952\n" +
                "M  CHG  1  11  -1\n" +
                "M  END";

        Chemical c = Chemical.parseMol(mol);
        IAtomContainer container = CdkUtil.toAtomContainer(c);
        sgroupsPropertyIsCollection(container);

        Chemical reparsed = Chemical.parseMol(c.toMol());
        IAtomContainer container2 = CdkUtil.toAtomContainer(reparsed);
        sgroupsPropertyIsCollection(container2);

        IAtomContainer clonedContainer = container.clone();
        sgroupsPropertyIsCollection(clonedContainer);
    }

    private void sgroupsPropertyIsCollection(IAtomContainer container) {
        Object groups = container.getProperty(CDKConstants.CTAB_SGROUPS);
        assertTrue(groups.getClass().toString(), groups instanceof Collection);
    }
}
