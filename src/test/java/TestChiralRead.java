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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import gov.nih.ncats.molwitch.cdk.CdkChemicalImpl;
import org.junit.Ignore;
import org.junit.Test;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.Chirality;
import gov.nih.ncats.molwitch.Stereocenter;
import gov.nih.ncats.molwitch.io.ChemicalReaderFactory;
import org.openscience.cdk.geometry.cip.CIPToolMod;

public class TestChiralRead {
		@Test
	   public void ensureChiralityIsTheSameRegardlessOfSDFOrMolOrAgnosticRead() throws Exception{
 	   
 	   String mm="\n" + 
					"  CDK     09262014552D\n" +
					"\n" + 
					" 69 71  0  0  1  0  0  0  0  0999 V2000\n" + 
					"    7.4100    3.6690    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.9100    3.6690    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.1600    2.3700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.9100    1.0710    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6600    2.3700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.9100    1.0700    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6600   -0.2200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.1600   -0.2175    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.9100   -1.5200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6602   -2.8189    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.1602   -2.8187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.9036   -1.5158    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.3784   -1.4938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.1792   -2.8143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    9.6791   -2.7943    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    7.4358   -4.1171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    5.8916   -4.0995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.4100   -1.5200    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6600   -2.8200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.4102   -4.1189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8400   -2.8200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5902   -4.1189    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5900   -1.5200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0900   -1.5200    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8400   -0.2200    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3400   -0.2200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.0900    1.0700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3400    2.3700    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8400    2.3700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0898    1.0711    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.0900    3.6700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -3.8402    4.9689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -5.3402    4.9687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.0904    6.2676    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -6.0900    3.6695    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5900    3.6700    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -0.8400    4.9700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -1.5902    6.2689    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6600    4.9700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.4102    6.2689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6604    7.5681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.4106    8.8670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    0.6608   10.1662    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.9106    8.8668    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    1.4100    3.6700    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    2.9100    3.6700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    3.6602    4.9689    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.5900    1.0675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -8.3422    2.3653    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -8.3378   -0.2328    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -7.7250   -1.5977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -8.8347   -2.6101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -10.1364   -1.8649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   -9.8338   -0.3903    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -10.8415    0.7208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -10.3831    2.1490    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -12.3076    0.4036    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -13.3153    1.5147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -12.8569    2.9429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -13.8647    4.0540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -15.3307    3.7369    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -13.4063    5.4823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -14.7814    1.1975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -15.2398   -0.2307    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -15.7891    2.3086    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -17.2552    1.9915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -18.2629    3.1025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -19.7290    2.7854    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  -17.8045    4.5308    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0\n" + 
					"  2  3  1  0  0  0  0\n" + 
					"  3  4  1  6  0  0  0\n" + 
					"  3  5  1  0  0  0  0\n" + 
					"  5  6  1  6  0  0  0\n" + 
					"  6  7  1  0  0  0  0\n" + 
					"  7  8  2  0  0  0  0\n" + 
					"  7  9  1  0  0  0  0\n" + 
					"  9 10  1  6  0  0  0\n" + 
					" 10 11  1  0  0  0  0\n" + 
					" 11 12  2  0  0  0  0\n" + 
					" 12 13  1  0  0  0  0\n" + 
					" 13 14  2  0  0  0  0\n" + 
					" 14 15  1  0  0  0  0\n" + 
					" 14 16  1  0  0  0  0\n" + 
					" 16 17  2  0  0  0  0\n" + 
					" 11 17  1  0  0  0  0\n" + 
					"  9 18  1  0  0  0  0\n" + 
					" 18 19  1  0  0  0  0\n" + 
					" 19 20  2  0  0  0  0\n" + 
					" 19 21  1  0  0  0  0\n" + 
					" 21 22  1  1  0  0  0\n" + 
					" 21 23  1  0  0  0  0\n" + 
					" 23 24  1  0  0  0  0\n" + 
					" 24 25  1  0  0  0  0\n" + 
					" 25 26  1  0  0  0  0\n" + 
					" 26 27  1  0  0  0  0\n" + 
					" 27 28  1  0  0  0  0\n" + 
					" 28 29  1  0  0  0  0\n" + 
					" 29 30  2  0  0  0  0\n" + 
					" 29 31  1  0  0  0  0\n" + 
					" 31 32  1  6  0  0  0\n" + 
					" 32 33  1  0  0  0  0\n" + 
					" 33 34  1  0  0  0  0\n" + 
					" 33 35  2  0  0  0  0\n" + 
					" 31 36  1  0  0  0  0\n" + 
					" 36 37  1  0  0  0  0\n" + 
					" 37 38  2  0  0  0  0\n" + 
					" 37 39  1  0  0  0  0\n" + 
					" 39 40  1  6  0  0  0\n" + 
					" 40 41  1  0  0  0  0\n" + 
					" 41 42  1  0  0  0  0\n" + 
					" 42 43  1  0  0  0  0\n" + 
					" 42 44  2  0  0  0  0\n" + 
					" 39 45  1  0  0  0  0\n" + 
					" 45 46  1  0  0  0  0\n" + 
					"  5 46  1  0  0  0  0\n" + 
					" 46 47  2  0  0  0  0\n" + 
					" 27 48  1  1  0  0  0\n" + 
					" 48 49  2  0  0  0  0\n" + 
					" 48 50  1  0  0  0  0\n" + 
					" 50 51  1  0  0  0  0\n" + 
					" 51 52  1  0  0  0  0\n" + 
					" 52 53  1  0  0  0  0\n" + 
					" 53 54  1  0  0  0  0\n" + 
					" 50 54  1  0  0  0  0\n" + 
					" 54 55  1  1  0  0  0\n" + 
					" 55 56  2  0  0  0  0\n" + 
					" 55 57  1  0  0  0  0\n" + 
					" 58 57  1  6  0  0  0\n" + 
					" 58 59  1  0  0  0  0\n" + 
					" 59 60  1  0  0  0  0\n" + 
					" 60 61  1  0  0  0  0\n" + 
					" 60 62  1  0  0  0  0\n" + 
					" 58 63  1  0  0  0  0\n" + 
					" 63 64  2  0  0  0  0\n" + 
					" 63 65  1  0  0  0  0\n" + 
					" 65 66  1  0  0  0  0\n" + 
					" 66 67  1  0  0  0  0\n" + 
					" 67 68  1  0  0  0  0\n" + 
					" 67 69  2  0  0  0  0\n" + 
					"M  END";
			Chemical c=Chemical.parse(mm);

			
			String genChiral = c.atoms()
			 .map(ca->ca.getChirality())
			 .filter(ch->!ch.equals(Chirality.Non_Chiral))
			 .map(ch->ch.toString())
			 .collect(Collectors.joining());

			
			c=Chemical.parseMol(mm);
			String molChiral = c.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());

			
			c=new Chemical(ChemicalReaderFactory.read("sdf", mm + "\n$$$$"));
			
			String sdfChiral = c.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());

			
			assertEquals(genChiral,molChiral);
			assertEquals(genChiral,sdfChiral);
			
		}
		
		@Test
	   	public void testAxialStereoMarkedR() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "  ChemDraw12062311172D\n"
					+ "\n"
					+ " 19 20  0  0  0  0  0  0  0  0999 V2000\n"
					+ "   -0.2855   -1.7339    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.4969   -1.7128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.8881   -2.3895    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.1525   -0.5075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.7189    0.2114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.9348   -0.4863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -0.6344    0.2326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.1102    0.9305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    2.3260    0.2326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.0679    0.9516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.0679   -0.4440    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.6766    1.6494    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.8925    0.9516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -0.6767    1.6706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.8926    0.9728    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -0.6767   -1.1207    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.8926   -0.4229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.1102    2.3895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -2.3260    0.2960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  2  0      \n"
					+ "  2  3  1  0      \n"
					+ "  2  4  1  0      \n"
					+ "  4  5  2  0      \n"
					+ "  4  6  1  0      \n"
					+ "  5  7  1  0      \n"
					+ "  5  8  1  6      \n"
					+ "  6  9  2  0      \n"
					+ "  7 10  2  0      \n"
					+ "  7 11  1  0      \n"
					+ "  8 12  1  0      \n"
					+ "  8 13  2  0      \n"
					+ "  9 13  1  0      \n"
					+ " 10 14  1  0      \n"
					+ " 10 15  1  0      \n"
					+ " 11 16  1  0      \n"
					+ " 11 17  2  0      \n"
					+ " 14 18  1  0      \n"
					+ " 15 19  2  0      \n"
					+ " 17 19  1  0      \n"
					+ "M  END");
			System.out.println(c1.toMol());
			Optional<Chirality> opChi=c1.atoms()
			  .filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.R, opChi.get());
	   	}
		@Test
	   	public void testAxialStereoMarkedS() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "  ChemDraw12062311172D\n"
					+ "\n"
					+ " 19 20  0  0  0  0  0  0  0  0999 V2000\n"
					+ "   -0.2855   -1.7339    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.4969   -1.7128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.8881   -2.3895    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.1525   -0.5075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.7189    0.2114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.9348   -0.4863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -0.6344    0.2326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.1102    0.9305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    2.3260    0.2326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.0679    0.9516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.0679   -0.4440    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    0.6766    1.6494    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "    1.8925    0.9516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -0.6767    1.6706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.8926    0.9728    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -0.6767   -1.1207    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.8926   -0.4229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -1.1102    2.3895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   -2.3260    0.2960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  2  0      \n"
					+ "  2  3  1  0      \n"
					+ "  2  4  1  0      \n"
					+ "  4  5  2  0      \n"
					+ "  4  6  1  0      \n"
					+ "  5  7  1  0      \n"
					+ "  5  8  1  1      \n"
					+ "  6  9  2  0      \n"
					+ "  7 10  2  0      \n"
					+ "  7 11  1  0      \n"
					+ "  8 12  1  0      \n"
					+ "  8 13  2  0      \n"
					+ "  9 13  1  0      \n"
					+ " 10 14  1  0      \n"
					+ " 10 15  1  0      \n"
					+ " 11 16  1  0      \n"
					+ " 11 17  2  0      \n"
					+ " 14 18  1  0      \n"
					+ " 15 19  2  0      \n"
					+ " 17 19  1  0      \n"
					+ "M  END");
			System.out.println(c1.toMol());
			Optional<Chirality> opChi=c1.atoms()
			  .filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.S, opChi.get());
	   	}
		
		@Ignore
		@Test
	   	public void testAxialStereoUndefinedMarkedAsCenter() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "   JSDraw212072312182D\n"
					+ "\n"
					+ "  6  6  0  0  0  0            999 V2000\n"
					+ "   17.4720   -8.8400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   16.1210   -8.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   16.1210   -6.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   18.8230   -8.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   18.8230   -6.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   17.4720   -5.7200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  2  0  0  0  0\n"
					+ "  2  3  1  0  0  0  0\n"
					+ "  1  4  1  0  0  0  0\n"
					+ "  4  5  2  0  0  0  0\n"
					+ "  5  6  1  0  0  0  0\n"
					+ "  6  3  2  0  0  0  0\n"
					+ "M  END");
			System.out.println(c1.toMol());
			Optional<Chirality> opChi=c1.atoms()
			  .filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.Parity_Either, opChi.get());
	   	}
		@Test
	   	public void testSimpleTetrahedralStereoMarked() throws Exception {
			Chemical c1=Chemical.parse("\n" + 
					"   JSDraw209282010112D\n" + 
					"\n" + 
					"  5  4  0  0  0  0              0 V2000\n" + 
					"   23.6688   -9.1798    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   25.1766   -9.2022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   26.5285   -9.1725    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   25.1766   -7.6422    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   25.1766  -10.7622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0\n" + 
					"  2  3  1  0  0  0  0\n" + 
					"  2  4  1  0  0  0  0\n" + 
					"  2  5  1  0  0  0  0\n" + 
					"M  END\n" + 
					"");
			
			Optional<Chirality> opChi=c1.atoms()
			  .filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.Parity_Either, opChi.get());
	   	}
		

		@Test
	   	public void testSimpleCisTransRingStereoMarked() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "   JSDraw212062312582D\n"
					+ "\n"
					+ "  8  8  0  0  1  0            999 V2000\n"
					+ "   24.4925   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   23.1415   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   23.1415   -5.7460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   24.4925   -4.9660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   25.8435   -5.7460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   25.8435   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   24.4925   -9.6460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   24.4925   -3.4060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  1  0  0  0  0\n"
					+ "  2  3  1  0  0  0  0\n"
					+ "  3  4  1  0  0  0  0\n"
					+ "  4  5  1  0  0  0  0\n"
					+ "  5  6  1  0  0  0  0\n"
					+ "  6  1  1  0  0  0  0\n"
					+ "  1  7  1  1  0  0  0\n"
					+ "  4  8  1  1  0  0  0\n"
					+ "M  END");

			assertEquals(2,c1.getAllStereocenters().size());
			Optional<Chirality> opChi=c1.atoms()
			  .filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.s, opChi.get());
	   	}
		
		@Test
	   	public void testSimpleCisTransRingStereoDetected() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "   JSDraw212062312582D\n"
					+ "\n"
					+ "  8  8  0  0  1  0            999 V2000\n"
					+ "   24.4925   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   23.1415   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   23.1415   -5.7460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   24.4925   -4.9660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   25.8435   -5.7460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   25.8435   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   24.4925   -9.6460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   24.4925   -3.4060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  1  0  0  0  0\n"
					+ "  2  3  1  0  0  0  0\n"
					+ "  3  4  1  0  0  0  0\n"
					+ "  4  5  1  0  0  0  0\n"
					+ "  5  6  1  0  0  0  0\n"
					+ "  6  1  1  0  0  0  0\n"
					+ "  1  7  1  0  0  0  0\n"
					+ "  4  8  1  0  0  0  0\n"
					+ "M  END");
			String sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("Parity_EitherParity_Either", sdfChiral);
			assertEquals(2,c1.getAllStereocenters().size());
//			String sdfChiral = c1.atoms()
//					 .map(ca->ca.getChirality())
//					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
//					 .map(ch->ch.toString())
//					 .collect(Collectors.joining());
//			
//			assertEquals("SSSR", sdfChiral);
	   	}
		
		
		@Test
	   	public void testSimpleTetrahedralStereoMarkedAsCenter() throws Exception {
			Chemical c1=Chemical.parse("\n" + 
					"   JSDraw209282010112D\n" + 
					"\n" + 
					"  5  4  0  0  0  0              0 V2000\n" + 
					"   23.6688   -9.1798    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   25.1766   -9.2022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   26.5285   -9.1725    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   25.1766   -7.6422    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   25.1766  -10.7622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0\n" + 
					"  2  3  1  0  0  0  0\n" + 
					"  2  4  1  0  0  0  0\n" + 
					"  2  5  1  0  0  0  0\n" + 
					"M  END\n" + 
					"");
			
			Optional<Chirality> opChi=c1.getAllStereocenters().stream()
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.Parity_Either, opChi.get());
	   	}
		
		@Test
	   	public void testSimpleTetrahedralStereoMarked2() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "   JSDraw212072312292D\n"
					+ "\n"
					+ " 42 42  0  0  1  0            999 V2000\n"
					+ "   33.4515  -13.1881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   32.1004  -13.9682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   30.7493  -13.1881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   39.1857  -26.0176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   40.4483  -25.1011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   40.5559  -23.5448    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   39.2082  -22.7591    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   41.9607  -22.8654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   42.1324  -21.3146    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   43.5022  -20.5673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   44.8574  -21.3399    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   43.3492  -19.0159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   42.1023  -18.0774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   41.7950  -16.5527    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   43.1529  -15.7973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   43.1828  -14.2274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   41.8398  -13.4353    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   41.9900  -11.8924    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   40.4817  -14.1926    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   40.4545  -15.7495    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   39.1482  -13.3918    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   39.2963  -11.8491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   37.9520  -11.0842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   37.8735   -9.5499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   36.5716   -8.7179    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   36.6213   -7.1692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   37.9730   -6.4054    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   39.0516   -7.5281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   40.5290   -8.0258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   42.0698   -7.7805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   38.8419   -2.1710    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   37.8550   -3.3778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   37.5422   -4.9049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   44.6029  -18.0888    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   44.9199  -16.5603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   44.6823  -15.0152    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   46.2997  -15.8321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   47.6046  -16.6791    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   46.3490  -14.2696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   47.5823  -13.3194    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   43.2388  -23.7596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   43.3766  -25.3138    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  1  0  0  0  0\n"
					+ "  2  3  1  0  0  0  0\n"
					+ "  4  5  1  0  0  0  0\n"
					+ "  5  6  1  0  0  0  0\n"
					+ "  6  7  1  6  0  0  0\n"
					+ "  6  8  1  0  0  0  0\n"
					+ "  8  9  1  6  0  0  0\n"
					+ "  9 10  1  0  0  0  0\n"
					+ " 10 11  2  0  0  0  0\n"
					+ " 10 12  1  0  0  0  0\n"
					+ " 12 13  1  0  0  0  0\n"
					+ " 13 14  1  0  0  0  0\n"
					+ " 14 15  1  0  0  0  0\n"
					+ " 15 16  1  0  0  0  0\n"
					+ " 16 17  1  0  0  0  0\n"
					+ " 17 18  2  0  0  0  0\n"
					+ " 17 19  1  0  0  0  0\n"
					+ " 19 20  1  0  0  0  0\n"
					+ " 14 20  1  0  0  0  0\n"
					+ " 19 21  1  0  0  0  0\n"
					+ " 21 22  2  0  0  0  0\n"
					+ " 22 23  1  0  0  0  0\n"
					+ " 23 24  1  0  0  0  0\n"
					+ " 24 25  1  0  0  0  0\n"
					+ " 25 26  2  0  0  0  0\n"
					+ " 26 27  1  0  0  0  0\n"
					+ " 27 28  1  0  0  0  0\n"
					+ " 28 29  1  0  0  0  0\n"
					+ " 29 30  1  0  0  0  0\n"
					+ " 32 31  1  0  0  0  0\n"
					+ " 32 33  1  0  0  0  0\n"
					+ " 27 33  1  0  0  0  0\n"
					+ " 12 34  1  6  0  0  0\n"
					+ " 34 35  1  0  0  0  0\n"
					+ " 35 36  2  0  0  0  0\n"
					+ " 35 37  1  0  0  0  0\n"
					+ " 37 38  1  6  0  0  0\n"
					+ " 37 39  1  0  0  0  0\n"
					+ " 39 40  1  0  0  0  0\n"
					+ "  8 41  1  0  0  0  0\n"
					+ " 41 42  2  0  0  0  0\n"
					+ " 31 30  2  0  0  0  0\n"
					+ "M  END");
			System.out.println(c1.toMol());
			String sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("SSSR", sdfChiral);
	   	}
		@Test
	   	public void testSetChiralityOnUnsetCaseWorksS() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "   JSDraw212062317442D\n"
					+ "\n"
					+ "  6  5  0  0  0  0            999 V2000\n"
					+ "   10.6080   -5.5640    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   11.9590   -4.7840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   13.3100   -5.5640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   14.6610   -4.7840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   16.0120   -5.5640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   11.9590   -3.2240    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  1  0  0  0  0\n"
					+ "  2  3  1  0  0  0  0\n"
					+ "  3  4  1  0  0  0  0\n"
					+ "  4  5  1  0  0  0  0\n"
					+ "  2  6  1  0  0  0  0\n"
					+ "M  END");
//			c1.generateCoordinates();
			String sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("Parity_Either", sdfChiral);
			c1.atoms()
			 .filter(ch->!ch.getChirality().equals(Chirality.Non_Chiral))
			 .forEach(ch->ch.setChirality(Chirality.S))
			 ;
			sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("S", sdfChiral);
			
			
	   	}
		@Test
	   	public void testSetChiralityOnUnsetCaseWorksR() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "   JSDraw212062317442D\n"
					+ "\n"
					+ "  6  5  0  0  0  0            999 V2000\n"
					+ "   10.6080   -5.5640    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   11.9590   -4.7840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   13.3100   -5.5640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   14.6610   -4.7840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   16.0120   -5.5640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   11.9590   -3.2240    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  1  0  0  0  0\n"
					+ "  2  3  1  0  0  0  0\n"
					+ "  3  4  1  0  0  0  0\n"
					+ "  4  5  1  0  0  0  0\n"
					+ "  2  6  1  0  0  0  0\n"
					+ "M  END");
//			c1.generateCoordinates();
			String sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("Parity_Either", sdfChiral);
			c1.atoms()
			 .filter(ch->!ch.getChirality().equals(Chirality.Non_Chiral))
			 .forEach(ch->ch.setChirality(Chirality.R))
			 ;
			sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("R", sdfChiral);
			
			
			
	   	}

		@Test
	   	public void testSulfoxideStereo() throws Exception {
			Chemical c1=Chemical.parse("\n" + 
					"   JSDraw209282010442D\n" + 
					"\n" + 
					"  5  4  0  0  1  0            999 V2000\n" + 
					"   16.7440   -8.8731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   18.0950   -8.0931    0.0000 S   0  3  0  0  0  0  0  0  0  0  0  0\n" + 
					"   19.4460   -8.8731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   20.7970   -8.0931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   18.0950   -6.5331    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0\n" + 
					"  2  3  1  0  0  0  0\n" + 
					"  3  4  1  0  0  0  0\n" + 
					"  2  5  1  1  0  0  0\n" + 
					"M  CHG  2   2   1   5  -1\n" + 
					"M  END");
			String sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("R", sdfChiral);
	   	}
		@Test
	   	public void testSulfoxideStereoPossible() throws Exception {
			Chemical c1=Chemical.parse("\n" + 
					"   JSDraw209282010462D\n" + 
					"\n" + 
					"  5  4  0  0  0  0            999 V2000\n" + 
					"   16.7440   -8.8731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   18.0950   -8.0931    0.0000 S   0  3  0  0  0  0  0  0  0  0  0  0\n" + 
					"   19.4460   -8.8731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   20.7970   -8.0931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   18.0950   -6.5331    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0\n" + 
					"  2  3  1  0  0  0  0\n" + 
					"  3  4  1  0  0  0  0\n" + 
					"  2  5  1  0  0  0  0\n" + 
					"M  CHG  2   2   1   5  -1\n" + 
					"M  END");
			String sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("Parity_Either", sdfChiral);
	   	}
		@Test
	   	public void testRemoveNonDescriptHydrogensDoesntRemoveStereoInformationOnMol() throws Exception {

	   		Chemical mol=Chemical.parse("\n" + 
	   				"   JSDraw209282016242D\n" + 
	   				"\n" + 
	   				"  8  7  0  0  1  0            999 V2000\n" + 
	   				"   28.5600   -9.6110    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   28.0990   -8.1210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   27.6380   -9.6110    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   29.4510   -7.3410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   30.9720   -7.6870    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   30.5110   -8.4850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   29.4510   -5.7810    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   26.7480   -7.3410    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"  1  2  1  0  0  0  0\n" + 
	   				"  2  3  1  1  0  0  0\n" + 
	   				"  2  4  1  0  0  0  0\n" + 
	   				"  4  5  1  6  0  0  0\n" + 
	   				"  4  6  1  0  0  0  0\n" + 
	   				"  4  7  1  0  0  0  0\n" + 
	   				"  2  8  1  0  0  0  0\n" + 
	   				"M  END");
	   		mol= Chemical.parse(mol.toMol());
	   		Chemical mol2= Chemical.createFromSmiles(mol.toSmiles());
	   		if(!mol2.hasCoordinates()){
	   			mol2.generateCoordinates();
	   		}

	   		mol2.removeNonDescriptHydrogens();
	   		mol2.generateCoordinates();
	   		mol2=Chemical.parse(mol2.toMol());
	   		
	   		String sdfChiral = mol2.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("SS", sdfChiral);
	   	}
		
		@Test
	   	public void testPsuedoStereocenterIsLowerCase() throws Exception {

	   		Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw212062315202D\n"
	   				+ "\n"
	   				+ " 10  9  0  0  1  0            999 V2000\n"
	   				+ "   14.4560   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   15.8070   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   15.8070   -5.7720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   17.1580   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   18.5090   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   14.4560   -9.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   13.1050   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   17.1580   -9.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   20.0690   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   11.7540   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  2  3  1  1  0  0  0\n"
	   				+ "  2  4  1  0  0  0  0\n"
	   				+ "  4  5  1  0  0  0  0\n"
	   				+ "  1  6  1  1  0  0  0\n"
	   				+ "  1  7  1  0  0  0  0\n"
	   				+ "  4  8  1  1  0  0  0\n"
	   				+ "  5  9  1  0  0  0  0\n"
	   				+ "  7 10  1  0  0  0  0\n"
	   				+ "M  END");
	   		mol= Chemical.parse(mol.toMol());
	   			   		assertEquals(3,mol.getAllStereocenters().size());
	   		String sdfChiral =mol.getAllStereocenters()
	   		.stream()
	   		.filter(Stereocenter::isDefined)
			.map(Stereocenter::getCenterAtom)
			.filter(a -> a.getChirality() != null)
			.map(a->a.getChirality())
			.map(ch->ch.toString())
			 .collect(Collectors.joining());
			
			assertEquals("SrR", sdfChiral);
	   	}
		

		@Test
	   	public void testPsuedoStereocenterIn135trimethylcyclohexaneIsLowerCase() throws Exception {

	   		Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw212082315262D\n"
	   				+ "\n"
	   				+ "  9  9  0  0  1  0            999 V2000\n"
	   				+ "   24.5440  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1930   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1930   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.5440   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.8950   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.8950   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.5440   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   21.8420  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   27.2460  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  2  3  1  0  0  0  0\n"
	   				+ "  3  4  1  0  0  0  0\n"
	   				+ "  4  5  1  0  0  0  0\n"
	   				+ "  5  6  1  0  0  0  0\n"
	   				+ "  6  1  1  0  0  0  0\n"
	   				+ "  4  7  1  1  0  0  0\n"
	   				+ "  2  8  1  1  0  0  0\n"
	   				+ "  6  9  1  6  0  0  0\n"
	   				+ "M  END");
	   		mol= Chemical.parse(mol.toMol());
	   			   		assertEquals(1,mol.getAllStereocenters().size());
	   		String sdfChiral =mol.getAllStereocenters()
	   		.stream()
	   		.filter(Stereocenter::isDefined)
			.map(Stereocenter::getCenterAtom)
			.filter(a -> a.getChirality() != null)
			.map(a->a.getChirality())
			.map(ch->ch.toString())
			 .collect(Collectors.joining());
			
			assertEquals("r", sdfChiral);
	   	}
		
		@Test
	   	public void testPsuedoStereocenterIn135trimethylcyclohexaneIsLowerCaseC() throws Exception {

			Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw212082315262D\n"
	   				+ "\n"
	   				+ "  9  9  0  0  1  0            999 V2000\n"
	   				+ "   24.5440  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1930   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1930   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.5440   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.8950   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.8950   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.5440   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   21.8420  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   27.2460  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  2  3  1  0  0  0  0\n"
	   				+ "  3  4  1  0  0  0  0\n"
	   				+ "  4  5  1  0  0  0  0\n"
	   				+ "  5  6  1  0  0  0  0\n"
	   				+ "  6  1  1  0  0  0  0\n"
	   				+ "  4  7  1  1  0  0  0\n"
	   				+ "  2  8  1  6  0  0  0\n"
	   				+ "  6  9  1  6  0  0  0\n"
	   				+ "M  END");
	   		mol= Chemical.parse(mol.toMol());
	   			   		assertEquals(1,mol.getAllStereocenters().size());
	   		String sdfChiral =mol.getAllStereocenters()
	   		.stream()
	   		.filter(Stereocenter::isDefined)
			.map(Stereocenter::getCenterAtom)
			.filter(a -> a.getChirality() != null)
			.map(a->a.getChirality())
			.map(ch->ch.toString())
			 .collect(Collectors.joining());
			
			assertEquals("r", sdfChiral);
	   	}
		
		@Test
	   	public void testPsuedoStereocenterIn135trimethylcyclohexaneIsLowerCaseB() throws Exception {

	   		Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw212082315262D\n"
	   				+ "\n"
	   				+ "  9  9  0  0  1  0            999 V2000\n"
	   				+ "   21.1640  -11.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   19.8130  -10.6600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   19.8130   -9.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   21.1640   -8.3200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   22.5150   -9.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   22.5150  -10.6600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   21.1640   -6.7600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   18.4620  -11.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.8660  -11.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  2  3  1  0  0  0  0\n"
	   				+ "  3  4  1  0  0  0  0\n"
	   				+ "  4  5  1  0  0  0  0\n"
	   				+ "  5  6  1  0  0  0  0\n"
	   				+ "  6  1  1  0  0  0  0\n"
	   				+ "  4  7  1  1  0  0  0\n"
	   				+ "  2  8  1  1  0  0  0\n"
	   				+ "  6  9  1  1  0  0  0\n"
	   				+ "M  END");
	   		mol= Chemical.parse(mol.toMol());
	   			   		assertEquals(3,mol.getAllStereocenters().size());
	   		String sdfChiral =mol.getAllStereocenters()
	   		.stream()
	   		.filter(Stereocenter::isDefined)
			.map(Stereocenter::getCenterAtom)
			.filter(a -> a.getChirality() != null)
			.map(a->a.getChirality())
			.map(ch->ch.toString())
			 .collect(Collectors.joining());
			
			assertEquals("sss", sdfChiral);
	   	}
		
		@Test
	   	public void testPsuedoStereocenterOnNonMesoNotFound() throws Exception {

	   		Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw212062315202D\n"
	   				+ "\n"
	   				+ " 10  9  0  0  1  0            999 V2000\n"
	   				+ "   14.4560   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   15.8070   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   15.8070   -5.7720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   17.1580   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   18.5090   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   14.4560   -9.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   13.1050   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   17.1580   -9.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   20.0690   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   11.7540   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  2  3  1  1  0  0  0\n"
	   				+ "  2  4  1  0  0  0  0\n"
	   				+ "  4  5  1  0  0  0  0\n"
	   				+ "  1  6  1  1  0  0  0\n"
	   				+ "  1  7  1  0  0  0  0\n"
	   				+ "  4  8  1  6  0  0  0\n"
	   				+ "  5  9  1  0  0  0  0\n"
	   				+ "  7 10  1  0  0  0  0\n"
	   				+ "M  END");
	   		mol= Chemical.parse(mol.toMol());
	   		assertEquals(2,mol.getAllStereocenters().size());
	   		String sdfChiral =mol.getAllStereocenters()
							   	 .stream()
							   	 .filter(Stereocenter::isDefined)
								 .map(Stereocenter::getCenterAtom)
								 .filter(a -> a.getChirality() != null)
								 .map(a->a.getChirality())
								 .map(ch->ch.toString())
								 .collect(Collectors.joining());
			
			assertEquals("SS", sdfChiral);
	   	}
		
		@Test
	   	public void testPsuedoStereocenterIsDetected() throws Exception {

	   		Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw212062315202D\n"
	   				+ "\n"
	   				+ " 10  9  0  0  1  0            999 V2000\n"
	   				+ "   14.4560   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   15.8070   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   15.8070   -5.7720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   17.1580   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   18.5090   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   14.4560   -9.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   13.1050   -7.3320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   17.1580   -9.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   20.0690   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   11.7540   -8.1120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  2  3  1  0  0  0  0\n"
	   				+ "  2  4  1  0  0  0  0\n"
	   				+ "  4  5  1  0  0  0  0\n"
	   				+ "  1  6  1  0  0  0  0\n"
	   				+ "  1  7  1  0  0  0  0\n"
	   				+ "  4  8  1  0  0  0  0\n"
	   				+ "  5  9  1  0  0  0  0\n"
	   				+ "  7 10  1  0  0  0  0\n"
	   				+ "M  END");
	   		mol= Chemical.parse(mol.toMol());
//	   		String sdfChiral = mol.atoms()
//					 .map(ca->ca.getChirality())
//					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
//					 .map(ch->ch.toString())
//					 .collect(Collectors.joining());
	   		assertEquals(3,mol.getAllStereocenters().size());
	   	}
		
		@Test
	   	public void testHavingExplicitHydrogenOnStereoCenterDoesNotInvalidatePhosphateCenter() throws Exception {

	   		Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw210312314162D\n"
	   				+ "\n"
	   				+ " 12 12  0  0  0  0            999 V2000\n"
	   				+ "   23.1415   -8.0860    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.4925   -8.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1415   -6.5260    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.8435   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1415   -9.6460    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.8435   -6.5260    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   27.1945   -5.7460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   28.5455   -6.5260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   28.5455   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   27.1945   -8.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   21.7905   -8.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.8435   -9.6460    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  1  3  2  0  0  0  0\n"
	   				+ "  2  4  1  0  0  0  0\n"
	   				+ "  1  5  1  0  0  0  0\n"
	   				+ "  4  6  1  0  0  0  0\n"
	   				+ "  6  7  1  0  0  0  0\n"
	   				+ "  7  8  1  0  0  0  0\n"
	   				+ "  8  9  1  0  0  0  0\n"
	   				+ "  9 10  1  0  0  0  0\n"
	   				+ " 10  4  1  0  0  0  0\n"
	   				+ "  1 11  1  0  0  0  0\n"
	   				+ "  4 12  1  0  0  0  0\n"
	   				+ "M  END");
	   		mol= Chemical.parse(mol.toMol());
	   		mol.removeNonDescriptHydrogens();
	   		System.out.println(mol.toMol());
	   		assertEquals(2,mol.getTetrahedrals().size());
	   	}
		
		@Test
	   	public void testHavingBondTableOrderChangedShouldKeepSameStereo() throws Exception {

			
	   		Chemical mol=Chemical.parse("\n"
	   				+ "   JSDraw212042310502D\n"
	   				+ "\n"
	   				+ " 14 14  0  0  1  0            999 V2000\n"
	   				+ "   21.4047   -8.8605    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.0654   -9.2811    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.4188   -8.5056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1278   -7.7212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.9094   -6.8325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   26.0050   -8.5767    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   26.6898   -5.7439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   28.1943   -6.0199    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.4240   -6.9456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.5922   -5.3944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   21.7756   -6.0051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   26.5338   -4.4447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   20.6247  -10.2116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.0654  -10.8412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  3  6  1  0  0  0  0\n"
	   				+ "  2  3  1  0  0  0  0\n"
	   				+ "  5  7  1  0  0  0  0\n"
	   				+ "  7  8  1  0  0  0  0\n"
	   				+ "  3  9  1  0  0  0  0\n"
	   				+ "  9  4  1  0  0  0  0\n"
	   				+ "  5  9  1  0  0  0  0\n"
	   				+ "  9 10  1  1  0  0  0\n"
	   				+ "  4 11  1  0  0  0  0\n"
	   				+ "  7 12  1  0  0  0  0\n"
	   				+ "  5  6  1  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  1 13  1  0  0  0  0\n"
	   				+ "  2 14  1  0  0  0  0\n"
	   				+ "M  END");
			Optional<Chirality> opChi=mol.getTetrahedrals().stream()
			  .map(ca->ca.getChirality())
			  .filter(ca->!ca.isEither())
			  			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.R, opChi.get());
			

	   		mol=Chemical.parse("\n"
	   				+ "   JSDraw212042310502D\n"
	   				+ "\n"
	   				+ " 14 14  0  0  1  0            999 V2000\n"
	   				+ "   21.4047   -8.8605    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.0654   -9.2811    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.4188   -8.5056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.1278   -7.7212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   25.9094   -6.8325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   26.0050   -8.5767    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   26.6898   -5.7439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   28.1943   -6.0199    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.4240   -6.9456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   24.5922   -5.3944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   21.7756   -6.0051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   26.5338   -4.4447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   20.6247  -10.2116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "   23.0654  -10.8412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
	   				+ "  2  3  1  0  0  0  0\n"
	   				+ "  3  6  1  0  0  0  0\n"
	   				+ "  5  7  1  0  0  0  0\n"
	   				+ "  7  8  1  0  0  0  0\n"
	   				+ "  3  9  1  0  0  0  0\n"
	   				+ "  9  4  1  0  0  0  0\n"
	   				+ "  5  9  1  0  0  0  0\n"
	   				+ "  9 10  1  1  0  0  0\n"
	   				+ "  4 11  1  0  0  0  0\n"
	   				+ "  7 12  1  0  0  0  0\n"
	   				+ "  5  6  1  0  0  0  0\n"
	   				+ "  1  2  1  0  0  0  0\n"
	   				+ "  1 13  1  0  0  0  0\n"
	   				+ "  2 14  1  0  0  0  0\n"
	   				+ "M  END");
			opChi=mol.getTetrahedrals().stream()
			  .map(ca->ca.getChirality())
			  .filter(ca->!ca.isEither())
			  			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.R, opChi.get());
	   	}
		
		@Test
	   	public void testReadingSmilesShouldNotGenerateCoordinates() throws Exception {

	   		Chemical mol=Chemical.parse("\n" + 
	   				"   JSDraw209282016242D\n" + 
	   				"\n" + 
	   				"  8  7  0  0  1  0            999 V2000\n" + 
	   				"   28.5600   -9.6110    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   28.0990   -8.1210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   27.6380   -9.6110    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   29.4510   -7.3410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   30.9720   -7.6870    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   30.5110   -8.4850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   29.4510   -5.7810    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   26.7480   -7.3410    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"  1  2  1  0  0  0  0\n" + 
	   				"  2  3  1  1  0  0  0\n" + 
	   				"  2  4  1  0  0  0  0\n" + 
	   				"  4  5  1  6  0  0  0\n" + 
	   				"  4  6  1  0  0  0  0\n" + 
	   				"  4  7  1  0  0  0  0\n" + 
	   				"  2  8  1  0  0  0  0\n" + 
	   				"M  END");
	   		Chemical mol2= Chemical.createFromSmiles(mol.toSmiles());
	   		
			
			assertFalse(mol2.hasCoordinates());
	   	}
		
	  	
	  	@Test
	   	public void testMethaneRemoveNonDescriptHydrogensMakesRightSmiles() throws Exception {
	   		Chemical mol=Chemical.parse("\n" + 
	   				"   JSDraw209282013552D\n" + 
	   				"\n" + 
	   				"  4  3  0  0  0  0            999 V2000\n" + 
	   				"   18.7200  -10.7328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   20.0710   -9.9528    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   18.7200  -12.2928    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"   17.3690   -9.9528    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
	   				"  1  2  1  0  0  0  0\n" + 
	   				"  1  3  1  0  0  0  0\n" + 
	   				"  1  4  1  0  0  0  0\n" + 
	   				"M  END");
	   		mol=Chemical.parse(mol.toSmiles());
	   		mol.removeNonDescriptHydrogens();
	   		assertEquals("C", mol.toSmiles());
	   	}
	  	
	  	@Test
	   	public void testSimpleExtendedTetrahedralStereoMarkedAsCenter() throws Exception {
			Chemical c1=Chemical.parse("\n" + 
					"   JSDraw209282018142D\n" + 
					"\n" + 
					"  7  6  0  0  1  0            999 V2000\n" + 
					"   10.0360   -7.4984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   11.3870   -6.7184    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   13.0510   -5.7824    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   10.0360   -9.0584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"    8.6850   -6.7184    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   13.0510   -4.2224    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   14.4540   -6.5624    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  2  0  0  0  0\n" + 
					"  2  3  2  0  0  0  0\n" + 
					"  1  4  1  0  0  0  0\n" + 
					"  1  5  1  0  0  0  0\n" + 
					"  3  6  1  6  0  0  0\n" + 
					"  3  7  1  1  0  0  0\n" + 
					"M  END");


			Optional<Chirality> opChi=c1.getAllStereocenters().stream()
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.S, opChi.get());
	   	}
	  	
		@Test
	   	public void testSimpleTetrahedralStereoOnQuatAmineMarkedAsCenter() throws Exception {
			Chemical c1=Chemical.parse("\n" + 
					"   JSDraw209282021202D\n" + 
					"\n" + 
					"  8  7  0  0  0  0              0 V2000\n" + 
					"   14.8200  -10.2960    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   16.1710   -9.5160    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n" + 
					"   16.1710   -7.9560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   17.5220  -10.2960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   14.8200   -7.1760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   18.8730   -9.5160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   18.8730   -7.9560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   16.1710  -11.0760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"  1  2  1  0  0  0  0\n" + 
					"  2  3  1  0  0  0  0\n" + 
					"  2  4  1  0  0  0  0\n" + 
					"  3  5  1  0  0  0  0\n" + 
					"  4  6  1  0  0  0  0\n" + 
					"  6  7  1  0  0  0  0\n" + 
					"  2  8  1  0  0  0  0\n" + 
					"M  CHG  1   2   1\n" + 
					"M  END\n" + 
					"");

			Optional<Chirality> opChi=c1.getAllStereocenters().stream()
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertTrue(opChi.isPresent());
			assertEquals(Chirality.Parity_Either, opChi.get());
	   	}
		
		@Test
	   	public void testSingleStereoBondBetweenTwoChiralCentersOnlyDefinesOneCenter() throws Exception {
			Chemical c1=Chemical.parse("\n"
					+ "   JSDraw210312313582D\n"
					+ "\n"
					+ "  9  9  0  0  1  0              0 V2000\n"
					+ "   17.6280  -11.3880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   16.2770  -10.6080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   16.2770   -9.0480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   17.6280   -8.2680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   18.9790   -9.0480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   18.9790  -10.6080    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   20.3300   -8.2680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   20.3300   -6.7080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "   21.6810   -9.0480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
					+ "  1  2  1  0  0  0  0\n"
					+ "  2  3  1  0  0  0  0\n"
					+ "  3  4  1  0  0  0  0\n"
					+ "  4  5  1  0  0  0  0\n"
					+ "  5  6  1  0  0  0  0\n"
					+ "  6  1  1  0  0  0  0\n"
					+ "  5  7  1  6  0  0  0\n"
					+ "  7  8  1  0  0  0  0\n"
					+ "  7  9  1  0  0  0  0\n"
					+ "M  END\n"
					+ "");
			
			assertEquals(2, c1.getAllStereocenters().size());
			assertEquals(Chirality.R, c1.getAllStereocenters().get(0).getChirality());
			assertEquals(Chirality.Parity_Either, c1.getAllStereocenters().get(1).getChirality());

	   	}
		
		@Test
	   	public void testQueryStructureSimpleTetrahedralStereoDoesntError() throws Exception {
			Chemical c1=Chemical.parse("S(=O)(=O)(O)O-C-[#6]");

			Optional<Chirality> opChi=c1.getAllStereocenters().stream()
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertFalse(opChi.isPresent());
	   	}

	@Test
	public void testSimpleChirality() throws Exception {
		Chemical c1=Chemical.parse("\n" +
				"  \n" +
				"\n" +
				"  5  4  0  0  1  0  0  0  0  0999 V2000\n" +
				"    6.4688   -5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    7.4916   -4.4094    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
				"    7.4916   -3.2280    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5148   -5.0002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.5380   -4.4094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"  1  2  1  0  0  0  0\n" +
				"  2  3  1  1  0  0  0\n" +
				"  2  4  1  0  0  0  0\n" +
				"  4  5  1  0  0  0  0\n" +
				"M  END\n");
		System.out.println(c1.toMol());
		Optional<Chirality> opChi=c1.atoms()
				.filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
				.map(ca->ca.getChirality())
				.findFirst();
		assertTrue(opChi.isPresent());
		assertEquals(Chirality.R, opChi.get());
	}

	@Test
	public void testSlowChiralityFalse() throws Exception {
		Chemical c1=Chemical.parse("\n" +
				"   JSDraw203122416262D\n" +
				"\n" +
				" 94104  0  0  0  0              0 V2000\n" +
				"   17.9069   -5.6149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   16.4061   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   23.3259   -5.5592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   21.7975   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   28.6894   -5.6149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   27.1888   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   32.5801   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   34.1364   -5.8788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730   -6.1554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730   -4.5813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5251   -6.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.4217   -5.4621    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.5439  -10.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   21.7975   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.0723  -10.1724    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   16.4061   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9937   -5.2968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   24.6875   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.0693   -5.4184    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.9048   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3711   -5.2462    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   30.0510   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6725   -4.5311    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   26.8171   -7.8935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   32.2085   -7.8935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   21.4269   -7.8935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   16.0429   -7.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   19.7188   -7.9213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   19.2685   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   30.5001   -7.9213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.9048   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   19.2685   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   30.0510   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   27.1888   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   32.5523   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9328  -10.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   24.6875   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3223  -10.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   18.0178  -10.0335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   23.4093  -10.0335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   28.8006  -10.0335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.2155   -7.7888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6258   -6.1554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6725  -10.8178    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6725   -9.2578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.4585   -9.7479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   34.0807   -9.7479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.1555   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3223   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730  -10.8178    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730   -9.2578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9328   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.5439   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.1555  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.5439  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9328  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3223  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.3209   -5.8788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   35.3095   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   35.3095   -6.9756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.2208   -6.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.0012   -6.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.2208   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.0012   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5251   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.2996   -7.7888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.0235  -11.5978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   39.3745  -10.8178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.7255  -11.5978    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.1025   -3.7593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   39.4858   -4.6118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.9158   -3.8401    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5098   -3.7942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    7.1466   -4.5813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.7834   -3.7942    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5220  -11.5978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    7.1710  -10.8178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.8200  -11.5978    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.9158   -2.2801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.9158   -5.4001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.7255  -10.0378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.7255  -13.1578    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.8200  -13.1578    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.8200  -10.0378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.7834   -5.3542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.7834   -2.2342    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    4.2234   -3.7942    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"    4.2600  -11.5978    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"   42.4758   -3.8401    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"   42.2855  -11.5978    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"  2 27  1  0  0  0  0\n" +
				" 58 61  1  0  0  0  0\n" +
				" 14 40  1  0  0  0  0\n" +
				" 15 16  1  0  0  0  0\n" +
				" 33 38  1  0  0  0  0\n" +
				" 33 41  1  0  0  0  0\n" +
				" 12 29  1  0  0  0  0\n" +
				" 61 63  1  0  0  0  0\n" +
				" 24 34  1  0  0  0  0\n" +
				" 47 59  1  0  0  0  0\n" +
				" 23 43  1  0  0  0  0\n" +
				" 36 56  2  0  0  0  0\n" +
				"  5 22  1  0  0  0  0\n" +
				" 50 51  1  0  0  0  0\n" +
				" 11 65  2  0  0  0  0\n" +
				" 15 31  1  0  0  0  0\n" +
				" 21 22  1  0  0  0  0\n" +
				" 31 46  1  0  0  0  0\n" +
				" 14 26  1  0  0  0  0\n" +
				"  3 18  1  0  0  0  0\n" +
				" 17 18  1  0  0  0  0\n" +
				" 26 28  1  0  0  0  0\n" +
				" 36 37  1  0  0  0  0\n" +
				" 13 55  2  0  0  0  0\n" +
				" 35 47  1  0  0  0  0\n" +
				"  9 61  2  0  0  0  0\n" +
				" 45 59  1  0  0  0  0\n" +
				" 28 32  1  0  0  0  0\n" +
				" 28 29  1  0  0  0  0\n" +
				"  3  4  1  0  0  0  0\n" +
				" 38 57  2  0  0  0  0\n" +
				"  7  8  1  0  0  0  0\n" +
				" 35 38  1  0  0  0  0\n" +
				" 59 60  2  0  0  0  0\n" +
				" 32 39  1  0  0  0  0\n" +
				"  6 17  1  0  0  0  0\n" +
				" 19 48  2  0  0  0  0\n" +
				" 62 64  1  0  0  0  0\n" +
				" 45 64  2  0  0  0  0\n" +
				" 25 35  1  0  0  0  0\n" +
				" 22 30  1  0  0  0  0\n" +
				" 13 32  1  0  0  0  0\n" +
				" 19 20  1  0  0  0  0\n" +
				" 16 27  1  0  0  0  0\n" +
				" 51 63  2  0  0  0  0\n" +
				" 15 54  2  0  0  0  0\n" +
				" 51 65  1  0  0  0  0\n" +
				" 21 49  2  0  0  0  0\n" +
				"  6 24  1  0  0  0  0\n" +
				" 18 42  1  0  0  0  0\n" +
				"  9 11  1  0  0  0  0\n" +
				" 30 33  1  0  0  0  0\n" +
				" 27 66  1  0  0  0  0\n" +
				" 43 60  1  0  0  0  0\n" +
				" 16 39  1  0  0  0  0\n" +
				"  4 12  1  0  0  0  0\n" +
				" 43 62  2  0  0  0  0\n" +
				" 24 42  1  0  0  0  0\n" +
				" 31 66  1  0  0  0  0\n" +
				"  4 26  1  0  0  0  0\n" +
				" 13 14  1  0  0  0  0\n" +
				"  8 60  1  0  0  0  0\n" +
				" 34 36  1  0  0  0  0\n" +
				" 17 52  2  0  0  0  0\n" +
				"  1 29  1  0  0  0  0\n" +
				" 34 41  1  0  0  0  0\n" +
				" 46 63  1  0  0  0  0\n" +
				" 25 30  1  0  0  0  0\n" +
				"  2 19  1  0  0  0  0\n" +
				"  1  2  1  0  0  0  0\n" +
				" 44 45  1  0  0  0  0\n" +
				"  5  6  1  0  0  0  0\n" +
				" 37 42  1  0  0  0  0\n" +
				" 12 53  2  0  0  0  0\n" +
				"  9 10  1  0  0  0  0\n" +
				" 37 40  1  0  0  0  0\n" +
				"  7 25  1  0  0  0  0\n" +
				" 20 58  1  0  0  0  0\n" +
				"  7 21  1  0  0  0  0\n" +
				" 20 66  1  0  0  0  0\n" +
				" 44 67  1  0  0  0  0\n" +
				" 67 68  1  0  0  0  0\n" +
				" 68 69  1  0  0  0  0\n" +
				" 23 70  1  0  0  0  0\n" +
				" 70 71  1  0  0  0  0\n" +
				" 71 72  1  0  0  0  0\n" +
				" 10 73  1  0  0  0  0\n" +
				" 73 74  1  0  0  0  0\n" +
				" 74 75  1  0  0  0  0\n" +
				" 50 76  1  0  0  0  0\n" +
				" 76 77  1  0  0  0  0\n" +
				" 77 78  1  0  0  0  0\n" +
				" 72 79  2  0  0  0  0\n" +
				" 72 80  2  0  0  0  0\n" +
				" 69 81  2  0  0  0  0\n" +
				" 69 82  2  0  0  0  0\n" +
				" 78 83  2  0  0  0  0\n" +
				" 78 84  2  0  0  0  0\n" +
				" 75 85  2  0  0  0  0\n" +
				" 75 86  2  0  0  0  0\n" +
				" 75 87  1  0  0  0  0\n" +
				" 78 88  1  0  0  0  0\n" +
				" 72 89  1  0  0  0  0\n" +
				" 69 90  1  0  0  0  0\n" +
				"M  STY  1   1 MUL\n" +
				"M  SAL   1  4  91  92  93  94\n" +
				"M  SPA   1  1  91\n" +
				"M  SDI   1  4   43.5776   -9.1000   43.5776   -6.7080\n" +
				"M  SDI   1  4   45.9176   -6.7080   45.9176   -9.1000\n" +
				"M  SMT   1 4\n" +
				"M  CHG  8  87  -1  88  -1  89  -1  90  -1  91   1  92   1  93   1  94   1\n" +
				"M  END\n");
		int ringCount = c1.getBondCount() - c1.getAtomCount() +1;
		System.out.printf("total atoms: %d bonds: %d; rings: %d\n", c1.getAtomCount(), c1.getBondCount(), ringCount);
		CdkChemicalImpl chem = (CdkChemicalImpl)c1.getImpl();

		CIPToolMod.label(chem.getContainer(), chem);
		List<Chirality> listChi=c1.atoms()
				.filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
				.map(ca->ca.getChirality())
				.collect(Collectors.toList());
		assertTrue(listChi.size() > 0);
		for (Chirality chirality : listChi) {
			System.out.printf("chirality: %s\n", chirality);
		}
	}

	@Test
	public void testSlowChiralityTrue() throws Exception {
		Chemical c1=Chemical.parse("\n" +
				"   JSDraw203122416262D\n" +
				"\n" +
				" 94104  0  0  0  0              0 V2000\n" +
				"   17.9069   -5.6149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   16.4061   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   23.3259   -5.5592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   21.7975   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   28.6894   -5.6149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   27.1888   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   32.5801   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   34.1364   -5.8788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730   -6.1554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730   -4.5813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5251   -6.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.4217   -5.4621    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.5439  -10.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   21.7975   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.0723  -10.1724    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   16.4061   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9937   -5.2968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   24.6875   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.0693   -5.4184    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.9048   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3711   -5.2462    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   30.0510   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6725   -4.5311    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   26.8171   -7.8935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   32.2085   -7.8935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   21.4269   -7.8935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   16.0429   -7.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   19.7188   -7.9213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   19.2685   -6.3095    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   30.5001   -7.9213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.9048   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   19.2685   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   30.0510   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   27.1888   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   32.5523   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9328  -10.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   24.6875   -9.2578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3223  -10.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   18.0178  -10.0335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   23.4093  -10.0335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   28.8006  -10.0335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.2155   -7.7888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6258   -6.1554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6725  -10.8178    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   36.6725   -9.2578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.4585   -9.7479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   34.0807   -9.7479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.1555   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3223   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730  -10.8178    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    9.8730   -9.2578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9328   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.5439   -3.8084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   15.1555  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   20.5439  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   25.9328  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   31.3223  -11.7565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.3209   -5.8788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   35.3095   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   35.3095   -6.9756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.2208   -6.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.0012   -6.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.2208   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.0012   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5251   -8.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.2996   -7.7888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.0235  -11.5978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   39.3745  -10.8178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.7255  -11.5978    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   38.1025   -3.7593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   39.4858   -4.6118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.9158   -3.8401    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5098   -3.7942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    7.1466   -4.5813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.7834   -3.7942    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    8.5220  -11.5978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    7.1710  -10.8178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.8200  -11.5978    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.9158   -2.2801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.9158   -5.4001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.7255  -10.0378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   40.7255  -13.1578    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.8200  -13.1578    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.8200  -10.0378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.7834   -5.3542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    5.7834   -2.2342    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"    4.2234   -3.7942    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"    4.2600  -11.5978    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"   42.4758   -3.8401    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"   42.2855  -11.5978    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"   44.3576   -7.9040    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n" +
				"  2 27  1  0  0  0  0\n" +
				" 58 61  1  0  0  0  0\n" +
				" 14 40  1  0  0  0  0\n" +
				" 15 16  1  0  0  0  0\n" +
				" 33 38  1  0  0  0  0\n" +
				" 33 41  1  0  0  0  0\n" +
				" 12 29  1  0  0  0  0\n" +
				" 61 63  1  0  0  0  0\n" +
				" 24 34  1  0  0  0  0\n" +
				" 47 59  1  0  0  0  0\n" +
				" 23 43  1  0  0  0  0\n" +
				" 36 56  2  0  0  0  0\n" +
				"  5 22  1  0  0  0  0\n" +
				" 50 51  1  0  0  0  0\n" +
				" 11 65  2  0  0  0  0\n" +
				" 15 31  1  0  0  0  0\n" +
				" 21 22  1  0  0  0  0\n" +
				" 31 46  1  0  0  0  0\n" +
				" 14 26  1  0  0  0  0\n" +
				"  3 18  1  0  0  0  0\n" +
				" 17 18  1  0  0  0  0\n" +
				" 26 28  1  0  0  0  0\n" +
				" 36 37  1  0  0  0  0\n" +
				" 13 55  2  0  0  0  0\n" +
				" 35 47  1  0  0  0  0\n" +
				"  9 61  2  0  0  0  0\n" +
				" 45 59  1  0  0  0  0\n" +
				" 28 32  1  0  0  0  0\n" +
				" 28 29  1  0  0  0  0\n" +
				"  3  4  1  0  0  0  0\n" +
				" 38 57  2  0  0  0  0\n" +
				"  7  8  1  0  0  0  0\n" +
				" 35 38  1  0  0  0  0\n" +
				" 59 60  2  0  0  0  0\n" +
				" 32 39  1  0  0  0  0\n" +
				"  6 17  1  0  0  0  0\n" +
				" 19 48  2  0  0  0  0\n" +
				" 62 64  1  0  0  0  0\n" +
				" 45 64  2  0  0  0  0\n" +
				" 25 35  1  0  0  0  0\n" +
				" 22 30  1  0  0  0  0\n" +
				" 13 32  1  0  0  0  0\n" +
				" 19 20  1  0  0  0  0\n" +
				" 16 27  1  0  0  0  0\n" +
				" 51 63  2  0  0  0  0\n" +
				" 15 54  2  0  0  0  0\n" +
				" 51 65  1  0  0  0  0\n" +
				" 21 49  2  0  0  0  0\n" +
				"  6 24  1  0  0  0  0\n" +
				" 18 42  1  0  0  0  0\n" +
				"  9 11  1  0  0  0  0\n" +
				" 30 33  1  0  0  0  0\n" +
				" 27 66  1  0  0  0  0\n" +
				" 43 60  1  0  0  0  0\n" +
				" 16 39  1  0  0  0  0\n" +
				"  4 12  1  0  0  0  0\n" +
				" 43 62  2  0  0  0  0\n" +
				" 24 42  1  0  0  0  0\n" +
				" 31 66  1  0  0  0  0\n" +
				"  4 26  1  0  0  0  0\n" +
				" 13 14  1  0  0  0  0\n" +
				"  8 60  1  0  0  0  0\n" +
				" 34 36  1  0  0  0  0\n" +
				" 17 52  2  0  0  0  0\n" +
				"  1 29  1  0  0  0  0\n" +
				" 34 41  1  0  0  0  0\n" +
				" 46 63  1  0  0  0  0\n" +
				" 25 30  1  0  0  0  0\n" +
				"  2 19  1  0  0  0  0\n" +
				"  1  2  1  0  0  0  0\n" +
				" 44 45  1  0  0  0  0\n" +
				"  5  6  1  0  0  0  0\n" +
				" 37 42  1  0  0  0  0\n" +
				" 12 53  2  0  0  0  0\n" +
				"  9 10  1  0  0  0  0\n" +
				" 37 40  1  0  0  0  0\n" +
				"  7 25  1  0  0  0  0\n" +
				" 20 58  1  0  0  0  0\n" +
				"  7 21  1  0  0  0  0\n" +
				" 20 66  1  0  0  0  0\n" +
				" 44 67  1  0  0  0  0\n" +
				" 67 68  1  0  0  0  0\n" +
				" 68 69  1  0  0  0  0\n" +
				" 23 70  1  0  0  0  0\n" +
				" 70 71  1  0  0  0  0\n" +
				" 71 72  1  0  0  0  0\n" +
				" 10 73  1  0  0  0  0\n" +
				" 73 74  1  0  0  0  0\n" +
				" 74 75  1  0  0  0  0\n" +
				" 50 76  1  0  0  0  0\n" +
				" 76 77  1  0  0  0  0\n" +
				" 77 78  1  0  0  0  0\n" +
				" 72 79  2  0  0  0  0\n" +
				" 72 80  2  0  0  0  0\n" +
				" 69 81  2  0  0  0  0\n" +
				" 69 82  2  0  0  0  0\n" +
				" 78 83  2  0  0  0  0\n" +
				" 78 84  2  0  0  0  0\n" +
				" 75 85  2  0  0  0  0\n" +
				" 75 86  2  0  0  0  0\n" +
				" 75 87  1  0  0  0  0\n" +
				" 78 88  1  0  0  0  0\n" +
				" 72 89  1  0  0  0  0\n" +
				" 69 90  1  0  0  0  0\n" +
				"M  STY  1   1 MUL\n" +
				"M  SAL   1  4  91  92  93  94\n" +
				"M  SPA   1  1  91\n" +
				"M  SDI   1  4   43.5776   -9.1000   43.5776   -6.7080\n" +
				"M  SDI   1  4   45.9176   -6.7080   45.9176   -9.1000\n" +
				"M  SMT   1 4\n" +
				"M  CHG  8  87  -1  88  -1  89  -1  90  -1  91   1  92   1  93   1  94   1\n" +
				"M  END\n");
		int ringCount = c1.getBondCount() - c1.getAtomCount() +1;
		System.out.printf("total atoms: %d bonds: %d; rings: %d\n", c1.getAtomCount(), c1.getBondCount(), ringCount);
		CdkChemicalImpl chem = (CdkChemicalImpl)c1.getImpl();

		CIPToolMod.label(chem.getContainer(), chem);
		List<Chirality> listChi=c1.atoms()
				.filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
				.map(ca->ca.getChirality())
				.collect(Collectors.toList());
		assertTrue(listChi.size() > 0);
		for (Chirality chirality : listChi) {
			System.out.printf("chirality: ");
		}
	}

	@Test
	public void testSlowChiralitySmall() throws Exception {
		String path = System.getenv(".");
		System.out.println("path: " + path);
		String userDirectory = new File("").getAbsolutePath();
		System.out.println("userDire: " + userDirectory);
		String fileName =userDirectory + "/src/test/resources/mols/pared_down.mol";
		File paredMol = new File(fileName);
		Chemical c1=Chemical.parse(Files.readString(paredMol.toPath()));
		int ringCount = c1.getBondCount() - c1.getAtomCount() +1;
		System.out.printf("total atoms: %d bonds: %d; rings: %d\n", c1.getAtomCount(), c1.getBondCount(), ringCount);
		CdkChemicalImpl chem = (CdkChemicalImpl)c1.getImpl();

		CIPToolMod.label(chem.getContainer(), chem);
		List<Chirality> listChi=c1.atoms()
				.filter(ca->ca.getChirality()!=Chirality.Non_Chiral)
				.map(ca->ca.getChirality())
				.collect(Collectors.toList());
		assertEquals(0, listChi.size());
		for (Chirality chirality : listChi) {
			System.out.printf("chirality: ");
		}
	}

	@Test
	public void testRingSystem() throws Exception {
		String path = System.getenv(".");
		System.out.println("path: " + path);
		String userDirectory = new File("").getAbsolutePath();
		System.out.println("userDirectory: " + userDirectory);
		String fileName =userDirectory + "/src/test/resources/mols/pared_down.mol";
		File paredMol = new File(fileName);
		CdkChemicalImpl c1= (CdkChemicalImpl) Chemical.parse(Files.readString(paredMol.toPath())).getImpl();
		CIPToolMod cipToolMod = new CIPToolMod();
		int ringSystemCount = cipToolMod.getSizeOfLargestRingSystem( c1);
		System.out.printf("total atoms: %d bonds: %d; rings: %d\n", c1.getAtomCount(), c1.getBondCount(), ringSystemCount);

		assertEquals(5, ringSystemCount);
	}

	@Test
	public void testRingSystem2() throws Exception {
		String path = System.getenv(".");
		System.out.println("path: " + path);
		String userDirectory = new File("").getAbsolutePath();
		System.out.println("userDirectory: " + userDirectory);
		String fileName =userDirectory + "/src/test/resources/mols/large.symmetric.mol";
		File paredMol = new File(fileName);
		CdkChemicalImpl c1= (CdkChemicalImpl) Chemical.parse(Files.readString(paredMol.toPath())).getImpl();
		CIPToolMod cipToolMod = new CIPToolMod();
		int ringSystemCount = cipToolMod.getSizeOfLargestRingSystem( c1);
		System.out.printf("total atoms: %d bonds: %d; rings: %d\n", c1.getAtomCount(), c1.getBondCount(), ringSystemCount);

		assertEquals(15, ringSystemCount);
	}

	@Test
	public void testRingSystem3() throws Exception {
		String path = System.getenv(".");
		System.out.println("path: " + path);
		String userDirectory = new File("").getAbsolutePath();
		System.out.println("userDirectory: " + userDirectory);
		String fileName =userDirectory + "/src/test/resources/mols/large.rings.pieces.mol";
		File paredMol = new File(fileName);
		CdkChemicalImpl c1= (CdkChemicalImpl) Chemical.parse(Files.readString(paredMol.toPath())).getImpl();
		CIPToolMod cipToolMod = new CIPToolMod();
		int ringSystemCount = cipToolMod.getSizeOfLargestRingSystem( c1);
		System.out.printf("total atoms: %d bonds: %d; rings: %d\n", c1.getAtomCount(), c1.getBondCount(), ringSystemCount);

		assertEquals(5, ringSystemCount);
	}

}
