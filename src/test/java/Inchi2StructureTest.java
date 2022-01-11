/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2022.
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
import gov.nih.ncats.molwitch.inchi.Inchi;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.layout.StructureDiagramGenerator;

import java.util.ArrayList;

public class Inchi2StructureTest {

    @Test
    public void inchiToStructure() throws Exception{
        String inchi = "InChI=1S/C16H25NO/c1-14(2)8-6-9-15(3)10-7-11-16(18)17-12-4-5-13-17/h7-8,10-11H,4-6,9,12-13H2,1-3H3/b11-7+,15-10+";
        Chemical c = Inchi.toChemical(inchi);

        System.out.println(c.getAtomCount());
    }

    @Test
    public void testGetInChIToStructure_String() throws CDKException {
        InChIToStructure parser = InChIGeneratorFactory.getInstance().getInChIToStructure(
                //     InChI=1S/C16H25NO/c1-14(2)8-6-9-15(3)10-7-11-16(18)17-12-4-5-13-17/h7-8,10-11H,4-6,9,12-13H2,1-3H3/b11-7+,15-10+
                "InChI=1S/C16H25NO/c1-14(2)8-6-9-15(3)10-7-11-16(18)17-12-4-5-13-17/h7-8,10-11H,4-6,9,12-13H2,1-3H3/b11-7+,15-10+",
                DefaultChemObjectBuilder.getInstance(), new ArrayList<>());

        StructureDiagramGenerator coordinateGenerator = new StructureDiagramGenerator(parser.getAtomContainer());
        coordinateGenerator.generateCoordinates();


    }
}
