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
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import java.util.Optional;
import java.util.stream.Collectors;

import org.junit.Test;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.Chirality;
import gov.nih.ncats.molwitch.TetrahedralChirality;
import gov.nih.ncats.molwitch.io.ChemicalReaderFactory;

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
			Chemical c1=Chemical.parse("\n" + 
					"   JSDraw209282010292D\n" + 
					"\n" + 
					" 26 26  0  0  1  0            999 V2000\n" + 
					"   41.5969   -8.3170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   40.0369   -8.3170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   39.2569   -9.6679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   40.0369  -11.0188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   37.6969   -9.6679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   36.9169  -11.0199    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   37.6969  -12.3614    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   39.2569  -12.3588    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   36.9169  -13.7134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   37.6971  -15.0642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   39.2571  -15.0640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   40.0302  -13.7090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   41.5640  -13.6862    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   42.3968  -15.0595    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   43.9567  -15.0387    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   41.6237  -16.4143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   40.0177  -16.3960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   35.3570  -13.7134    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   34.5770  -15.0654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   35.3572  -16.4162    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   33.0170  -15.0654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   32.2368  -16.4162    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   32.2370  -13.7134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   30.6771  -13.7134    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   36.9169   -8.3159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
					"   37.6971   -6.9651    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
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
					"  5 25  1  0  0  0  0\n" + 
					" 25 26  2  0  0  0  0\n" + 
					"M  END");
			String sdfChiral = c1.atoms()
					 .map(ca->ca.getChirality())
					 .filter(ch->!ch.equals(Chirality.Non_Chiral))
					 .map(ch->ch.toString())
					 .collect(Collectors.joining());
			
			assertEquals("SSSR", sdfChiral);
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
	   	public void testQueryStructureSimpleTetrahedralStereoDoesntError() throws Exception {
			Chemical c1=Chemical.parse("S(=O)(=O)(O)O-C-[#6]");

			Optional<Chirality> opChi=c1.getAllStereocenters().stream()
			  .map(ca->ca.getChirality())
			  .findFirst();
			assertFalse(opChi.isPresent());
	   	}
	  	
}
