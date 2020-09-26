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

import static org.junit.Assert.assertEquals;

import java.util.stream.Collectors;

import org.junit.Test;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.Chirality;
import gov.nih.ncats.molwitch.io.ChemicalReaderFactory;

public class TestChiralRead {
		@Test
	   public void ensureChiralityIsTheSameRegardlessOfSDFOrMolOrAgnosticRead() throws Exception{
 	   
 	   String mm="\n" + 
					"  CDK     09262014553D\n" + 
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
			
			System.out.println("CHIRALITY FROM GENERIC PARSE:" + genChiral);
			
			c=Chemical.parseMol(mm);
			String molChiral = c.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());

			System.out.println("CHIRALITY FROM MOL PARSE:" + molChiral);
			
			c=new Chemical(ChemicalReaderFactory.read("sdf", mm + "\n$$$$"));
			
			String sdfChiral = c.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			System.out.println("CHIRALITY FROM SDF PARSE:" + sdfChiral);
			
			assertEquals(genChiral,molChiral);
			assertEquals(genChiral,sdfChiral);
			
		}
}
