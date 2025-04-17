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

import static org.junit.Assert.assertEquals;

import java.util.stream.Collectors;

import org.junit.Test;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.Chirality;
import gov.nih.ncats.molwitch.Bond.BondType;
import gov.nih.ncats.molwitch.inchi.InChiResult;
import gov.nih.ncats.molwitch.io.ChemicalReaderFactory;

public class TestChemicalValidation {

	@Test
   	public void testPentavalent() throws Exception {
   		Chemical c=Chemical.parse("\n" + 
   				"   JSDraw209262021142D\n" + 
   				"\n" + 
   				"  6  5  0  0  0  0            999 V2000\n" + 
   				"   21.2160   -6.7742    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   22.5670   -5.9942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   21.2160   -8.3342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   19.8650   -5.9942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   21.2160   -5.2142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   19.8650   -7.5542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"  1  2  1  0  0  0  0\n" + 
   				"  1  3  1  0  0  0  0\n" + 
   				"  1  4  1  0  0  0  0\n" + 
   				"  1  5  1  0  0  0  0\n" + 
   				"  1  6  1  0  0  0  0\n" + 
   				"M  END");
   		
   		long terror=c.atoms()
   				.filter(ca->ca.hasValenceError()).count();
   		
   		assertEquals(1,terror);
   		
   		

   	}
 	@Test
   	public void testInchiAromatic() throws Exception {
   		Chemical c=Chemical.parse("\n" + 
   				"   JSDraw209262021582D\n" + 
   				"\n" + 
   				"  6  6  0  0  0  0              0 V2000\n" + 
   				"   21.8408   -6.9164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   23.1920   -6.1369    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   24.5432   -6.9164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   24.5432   -8.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   21.8408   -8.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   23.1920   -9.2551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"  1  2  4  0  0  0  0\n" + 
   				"  2  3  4  0  0  0  0\n" + 
   				"  3  4  4  0  0  0  0\n" + 
   				"  1  5  4  0  0  0  0\n" + 
   				"  5  6  4  0  0  0  0\n" + 
   				"  4  6  4  0  0  0  0\n" + 
   				"M  END");

        InChiResult result = c.toInchi();
        String key = result.getKey();
        assertEquals("UHOVQNZJYSORNB-UHFFFAOYSA-N", key);


   	}
   	
   	@Test
   	public void testSimplestDearomatization() throws Exception {
   		Chemical c=Chemical.parse("\n" + 
   				"   JSDraw209262021582D\n" + 
   				"\n" + 
   				"  6  6  0  0  0  0              0 V2000\n" + 
   				"   21.8408   -6.9164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   23.1920   -6.1369    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   24.5432   -6.9164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   24.5432   -8.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   21.8408   -8.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"   23.1920   -9.2551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
   				"  1  2  4  0  0  0  0\n" + 
   				"  2  3  4  0  0  0  0\n" + 
   				"  3  4  4  0  0  0  0\n" + 
   				"  1  5  4  0  0  0  0\n" + 
   				"  5  6  4  0  0  0  0\n" + 
   				"  4  6  4  0  0  0  0\n" + 
   				"M  END");
   		c.kekulize();
   		
   		long l=c.bonds().filter(b->b.getBondType().equals(BondType.SINGLE) || b.getBondType().equals(BondType.DOUBLE))
   		         .count();

   		assertEquals(6,l);


   	}
}
