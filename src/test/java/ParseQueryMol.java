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

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Optional;

import org.junit.Test;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.fingerprint.Fingerprint;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.FingerprintSpecification;
import gov.nih.ncats.molwitch.search.MolSearcherFactory;

public class ParseQueryMol {

    @Test
    public void atomListsAndQueryBondFingerprintsWork() throws IOException {
    	  String mol= "\n" + 
          		"   JSDraw209242017322D\n" + 
          		"\n" + 
          		" 22 22  0  0  0  0            999 V2000\n" + 
          		"   25.4825   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   24.1315   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   24.1315   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   26.8335   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   26.8335   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   25.4825   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   28.1845   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   28.1845   -3.0160    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   29.5355   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   22.7805   -4.5760    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   25.4825   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   26.8335  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   24.1315  -10.0360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   28.1845   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   29.5355  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   30.8865   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   32.2375  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   33.5885   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   32.2375  -11.5960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   33.5885   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   34.9395   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   33.5885  -12.3760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"  1  2  2  0  0  0  0\n" + 
          		"  2  3  1  0  0  0  0\n" + 
          		"  1  4  1  0  0  0  0\n" + 
          		"  4  5  2  0  0  0  0\n" + 
          		"  5  6  1  0  0  0  0\n" + 
          		"  6  3  2  0  0  0  0\n" + 
          		"  5  7  1  0  0  0  0\n" + 
          		"  7  8  1  0  0  0  0\n" + 
          		"  7  9  6  0  0  0  0\n" + 
          		"  3 10  1  0  0  0  0\n" + 
          		"  1 11  1  0  0  0  0\n" + 
          		" 11 12  1  0  0  0  0\n" + 
          		" 11 13  8  0  0  0  0\n" + 
          		" 12 14  1  0  0  0  0\n" + 
          		" 14 15  1  0  0  0  0\n" + 
          		" 15 16  1  0  0  0  0\n" + 
          		" 16 17  1  0  0  0  0\n" + 
          		" 17 18  1  0  0  0  0\n" + 
          		" 17 19  1  0  0  0  0\n" + 
          		" 18 20  1  0  0  0  0\n" + 
          		" 20 21  4  0  0  0  0\n" + 
          		" 19 22  7  0  0  0  0\n" + 
          		"M  ALS   8  2 F O   N   \n" + 
          		"M  ALS  10  1 T C   \n" + 
          		"M  END";
          
        
        Chemical c = Chemical.parseMol(mol);

        Fingerprinter fingerPrinterSub =  Fingerprinters.getFingerprinter(FingerprintSpecification.PATH_BASED.create().setLength(512));
        
        Fingerprint fp=fingerPrinterSub.computeFingerprint(c);
        
        
        String mm="\n" + 
        		"   JSDraw209242017082D\n" + 
        		"\n" + 
        		" 27 28  0  0  0  0            999 V2000\n" + 
        		"   25.4825   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -3.0160    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   22.7805   -4.5760    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315  -10.0360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   30.8865   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -11.5960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   34.9395   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885  -12.3760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335  -11.5960    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   34.9395   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"  1  2  2  0  0  0  0\n" + 
        		"  2  3  1  0  0  0  0\n" + 
        		"  1  4  1  0  0  0  0\n" + 
        		"  4  5  2  0  0  0  0\n" + 
        		"  5  6  1  0  0  0  0\n" + 
        		"  6  3  2  0  0  0  0\n" + 
        		"  5  7  1  0  0  0  0\n" + 
        		"  7  8  1  0  0  0  0\n" + 
        		"  7  9  1  0  0  0  0\n" + 
        		"  3 10  1  0  0  0  0\n" + 
        		"  1 11  1  0  0  0  0\n" + 
        		" 11 12  1  0  0  0  0\n" + 
        		" 11 13  1  0  0  0  0\n" + 
        		" 12 14  1  0  0  0  0\n" + 
        		" 14 15  1  0  0  0  0\n" + 
        		" 15 16  1  0  0  0  0\n" + 
        		" 16 17  1  0  0  0  0\n" + 
        		" 17 18  1  0  0  0  0\n" + 
        		" 17 19  1  0  0  0  0\n" + 
        		" 18 20  1  0  0  0  0\n" + 
        		" 20 21  4  0  0  0  0\n" + 
        		" 19 22  2  0  0  0  0\n" + 
        		" 12 23  1  0  0  0  0\n" + 
        		" 21 24  4  0  0  0  0\n" + 
        		" 20 25  4  0  0  0  0\n" + 
        		" 25 26  4  0  0  0  0\n" + 
        		" 26 27  4  0  0  0  0\n" + 
        		" 27 24  4  0  0  0  0\n" + 
        		"M  END";
        Chemical c2 = Chemical.parseMol(mm);

        Fingerprint fp2=fingerPrinterSub.computeFingerprint(c2);
        
        BitSet bsTemp = BitSet.valueOf(fp.toBitSet().toLongArray());
        bsTemp.and(fp2.toBitSet());

        assertEquals(fp.populationCount(), bsTemp.cardinality());
        
    }
    
    @Test
    public void atomListsAndQueryBondSmartsExportWorks() throws IOException {
    	  String mol= "\n" + 
          		"   JSDraw209242017322D\n" + 
          		"\n" + 
          		" 22 22  0  0  0  0            999 V2000\n" + 
          		"   25.4825   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   24.1315   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   24.1315   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   26.8335   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   26.8335   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   25.4825   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   28.1845   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   28.1845   -3.0160    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   29.5355   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   22.7805   -4.5760    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   25.4825   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   26.8335  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   24.1315  -10.0360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   28.1845   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   29.5355  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   30.8865   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   32.2375  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   33.5885   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   32.2375  -11.5960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   33.5885   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   34.9395   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"   33.5885  -12.3760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
          		"  1  2  2  0  0  0  0\n" + 
          		"  2  3  1  0  0  0  0\n" + 
          		"  1  4  1  0  0  0  0\n" + 
          		"  4  5  2  0  0  0  0\n" + 
          		"  5  6  1  0  0  0  0\n" + 
          		"  6  3  2  0  0  0  0\n" + 
          		"  5  7  1  0  0  0  0\n" + 
          		"  7  8  1  0  0  0  0\n" + 
          		"  7  9  6  0  0  0  0\n" + 
          		"  3 10  1  0  0  0  0\n" + 
          		"  1 11  1  0  0  0  0\n" + 
          		" 11 12  1  0  0  0  0\n" + 
          		" 11 13  8  0  0  0  0\n" + 
          		" 12 14  1  0  0  0  0\n" + 
          		" 14 15  1  0  0  0  0\n" + 
          		" 15 16  1  0  0  0  0\n" + 
          		" 16 17  1  0  0  0  0\n" + 
          		" 17 18  1  0  0  0  0\n" + 
          		" 17 19  1  0  0  0  0\n" + 
          		" 18 20  1  0  0  0  0\n" + 
          		" 20 21  4  0  0  0  0\n" + 
          		" 19 22  7  0  0  0  0\n" + 
          		"M  ALS   8  2 F O   N   \n" + 
          		"M  ALS  10  1 T C   \n" + 
          		"M  END";
          
        
        Chemical c = Chemical.parseMol(mol);
        String smarts=c.toSmarts();
        
        System.out.println(smarts);
        
        
        
        String mm="\n" + 
        		"   JSDraw209242017082D\n" + 
        		"\n" + 
        		" 27 28  0  0  0  0            999 V2000\n" + 
        		"   25.4825   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -3.0160    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   22.7805   -4.5760    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315  -10.0360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   30.8865   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -11.5960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   34.9395   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885  -12.3760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335  -11.5960    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   34.9395   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"  1  2  2  0  0  0  0\n" + 
        		"  2  3  1  0  0  0  0\n" + 
        		"  1  4  1  0  0  0  0\n" + 
        		"  4  5  2  0  0  0  0\n" + 
        		"  5  6  1  0  0  0  0\n" + 
        		"  6  3  2  0  0  0  0\n" + 
        		"  5  7  1  0  0  0  0\n" + 
        		"  7  8  1  0  0  0  0\n" + 
        		"  7  9  1  0  0  0  0\n" + 
        		"  3 10  1  0  0  0  0\n" + 
        		"  1 11  1  0  0  0  0\n" + 
        		" 11 12  1  0  0  0  0\n" + 
        		" 11 13  1  0  0  0  0\n" + 
        		" 12 14  1  0  0  0  0\n" + 
        		" 14 15  1  0  0  0  0\n" + 
        		" 15 16  1  0  0  0  0\n" + 
        		" 16 17  1  0  0  0  0\n" + 
        		" 17 18  1  0  0  0  0\n" + 
        		" 17 19  1  0  0  0  0\n" + 
        		" 18 20  1  0  0  0  0\n" + 
        		" 20 21  4  0  0  0  0\n" + 
        		" 19 22  2  0  0  0  0\n" + 
        		" 12 23  1  0  0  0  0\n" + 
        		" 21 24  4  0  0  0  0\n" + 
        		" 20 25  4  0  0  0  0\n" + 
        		" 25 26  4  0  0  0  0\n" + 
        		" 26 27  4  0  0  0  0\n" + 
        		" 27 24  4  0  0  0  0\n" + 
        		"M  END";
        Chemical c2 = Chemical.parseMol(mm);
        System.out.println(c2.toSmarts());

        
    }
   
    @Test
    public void atomListsAndQueryBondStructureSearchWork() throws IOException {
        String mol= "\n" + 
        		"   JSDraw209242017322D\n" + 
        		"\n" + 
        		" 22 22  0  0  0  0            999 V2000\n" + 
        		"   25.4825   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -3.0160    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   22.7805   -4.5760    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315  -10.0360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   30.8865   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -11.5960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   34.9395   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885  -12.3760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"  1  2  2  0  0  0  0\n" + 
        		"  2  3  1  0  0  0  0\n" + 
        		"  1  4  1  0  0  0  0\n" + 
        		"  4  5  2  0  0  0  0\n" + 
        		"  5  6  1  0  0  0  0\n" + 
        		"  6  3  2  0  0  0  0\n" + 
        		"  5  7  1  0  0  0  0\n" + 
        		"  7  8  1  0  0  0  0\n" + 
        		"  7  9  6  0  0  0  0\n" + 
        		"  3 10  1  0  0  0  0\n" + 
        		"  1 11  1  0  0  0  0\n" + 
        		" 11 12  1  0  0  0  0\n" + 
        		" 11 13  8  0  0  0  0\n" + 
        		" 12 14  1  0  0  0  0\n" + 
        		" 14 15  1  0  0  0  0\n" + 
        		" 15 16  1  0  0  0  0\n" + 
        		" 16 17  1  0  0  0  0\n" + 
        		" 17 18  1  0  0  0  0\n" + 
        		" 17 19  1  0  0  0  0\n" + 
        		" 18 20  1  0  0  0  0\n" + 
        		" 20 21  4  0  0  0  0\n" + 
        		" 19 22  7  0  0  0  0\n" + 
        		"M  ALS   8  2 F O   N   \n" + 
        		"M  ALS  10  1 T C   \n" + 
        		"M  END";
        
        Chemical c = Chemical.parseMol(mol);
        
        String mm="\n" + 
        		"   JSDraw209242017082D\n" + 
        		"\n" + 
        		" 27 28  0  0  0  0            999 V2000\n" + 
        		"   25.4825   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -3.0160    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   22.7805   -4.5760    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   25.4825   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   24.1315  -10.0360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   28.1845   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   29.5355  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   30.8865   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -10.0360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -9.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375  -11.5960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -7.6960    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   34.9395   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885  -12.3760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   26.8335  -11.5960    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   34.9395   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375   -6.9160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   32.2375   -5.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"   33.5885   -4.5760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
        		"  1  2  2  0  0  0  0\n" + 
        		"  2  3  1  0  0  0  0\n" + 
        		"  1  4  1  0  0  0  0\n" + 
        		"  4  5  2  0  0  0  0\n" + 
        		"  5  6  1  0  0  0  0\n" + 
        		"  6  3  2  0  0  0  0\n" + 
        		"  5  7  1  0  0  0  0\n" + 
        		"  7  8  1  0  0  0  0\n" + 
        		"  7  9  1  0  0  0  0\n" + 
        		"  3 10  1  0  0  0  0\n" + 
        		"  1 11  1  0  0  0  0\n" + 
        		" 11 12  1  0  0  0  0\n" + 
        		" 11 13  1  0  0  0  0\n" + 
        		" 12 14  1  0  0  0  0\n" + 
        		" 14 15  1  0  0  0  0\n" + 
        		" 15 16  1  0  0  0  0\n" + 
        		" 16 17  1  0  0  0  0\n" + 
        		" 17 18  1  0  0  0  0\n" + 
        		" 17 19  1  0  0  0  0\n" + 
        		" 18 20  1  0  0  0  0\n" + 
        		" 20 21  4  0  0  0  0\n" + 
        		" 19 22  2  0  0  0  0\n" + 
        		" 12 23  1  0  0  0  0\n" + 
        		" 21 24  4  0  0  0  0\n" + 
        		" 20 25  4  0  0  0  0\n" + 
        		" 25 26  4  0  0  0  0\n" + 
        		" 26 27  4  0  0  0  0\n" + 
        		" 27 24  4  0  0  0  0\n" + 
        		"M  END";
        
        System.out.println(mm);
        
        Chemical c2 = Chemical.parseMol(mm);
        
        Optional<int[]> hit = MolSearcherFactory.create(c).search(c2);
		

		assertTrue(hit.isPresent());
        
        
    }
    
       @Test
       public void ensureAnyBondAndAromatizationWorksForSubstructureSearch() throws IOException {
           String mol= "\n" + 
           		"   JSDraw209252000242D\n" + 
           		"\n" + 
           		"  8  8  0  0  0  0            999 V2000\n" + 
           		"   28.1035   -4.9661    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   29.4545   -5.7461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   30.8055   -4.9661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   32.1565   -5.7461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   32.1565   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   30.8055   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   29.4545   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   33.5073   -8.0860    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"  1  2  8  0  0  0  0\n" + 
           		"  2  3  1  0  0  0  0\n" + 
           		"  3  4  2  0  0  0  0\n" + 
           		"  4  5  1  0  0  0  0\n" + 
           		"  5  6  2  0  0  0  0\n" + 
           		"  6  7  1  0  0  0  0\n" + 
           		"  2  7  2  0  0  0  0\n" + 
           		"  5  8  1  0  0  0  0\n" + 
           		"M  ALS   1  2 F N   O   \n" + 
           		"M  END";
           
           
           Chemical c = Chemical.parse(mol);


           try{
        	   c.generateCoordinates();
           }catch(Exception e){
        	   e.printStackTrace();
           }
           c.aromatize();
           System.out.println(c.toMol());
           System.out.println(c.toSmarts());
           
           
           String mol2="NC1=CC(O)=C(O)C=C1";
           Chemical c2 = Chemical.parse(mol2);
           c2.aromatize();
           try{
        	   c2.generateCoordinates();
           }catch(Exception e){
        	   e.printStackTrace();
           }
           System.out.println(c2.toMol());
           System.out.println(c2.toSmiles());
           Optional<int[]> hit = MolSearcherFactory.create(c).search(c2);
           assertTrue(hit.isPresent());
       }
       
       
       @Test
       public void ensureAnyBondAndAromatizationInComplexExampleWorksForSubstructureSearch() throws IOException {
           String mol= "\n" + 
           		"  CDK     09252009303D\n" + 
           		"\n" + 
           		" 19 21  0  0  0  0  0  0  0  0999 V2000\n" + 
           		"   -1.9500    2.4000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -1.9500    0.9000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -0.6500    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    0.6500    0.9000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    1.9500    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    3.0600    1.1500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    2.4500    2.5200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    3.8762    2.9846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    2.1359    3.9867    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    0.9600    2.3600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -0.0451    3.4735    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    1.9500   -1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    0.6500   -2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"    0.6500   -3.6000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -0.6500   -1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -1.9500   -2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -1.9500   -3.6000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -3.2500   -1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   -3.2500    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"  1  2  8  0  0  0  0\n" + 
           		"  2 19  4  0  0  0  0\n" + 
           		"  2  3  4  0  0  0  0\n" + 
           		"  3 15  4  0  0  0  0\n" + 
           		"  3  4  4  0  0  0  0\n" + 
           		"  4 10  6  0  0  0  0\n" + 
           		"  4  5  4  0  0  0  0\n" + 
           		"  5  6  6  0  0  0  0\n" + 
           		"  6  7  6  0  0  0  0\n" + 
           		"  7  8  6  0  0  0  0\n" + 
           		"  7  9  6  0  0  0  0\n" + 
           		"  7 10  6  0  0  0  0\n" + 
           		" 10 11  2  0  0  0  0\n" + 
           		"  5 12  4  0  0  0  0\n" + 
           		" 12 13  4  0  0  0  0\n" + 
           		" 13 14  6  0  0  0  0\n" + 
           		" 13 15  4  0  0  0  0\n" + 
           		" 15 16  4  0  0  0  0\n" + 
           		" 16 17  6  0  0  0  0\n" + 
           		" 16 18  4  0  0  0  0\n" + 
           		" 18 19  4  0  0  0  0\n" + 
      		    "M  ALS   1  2 F O   N   \n" + 
           		"M  END";
                    
           Chemical c = Chemical.parse(mol);
           try{
        	   c.generateCoordinates();
           }catch(Exception e){
        	   e.printStackTrace();
           }
           c.aromatize();
           System.out.println(c.toMol());
           System.out.println(c.toSmarts());
           
           
           String mol2="COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12";
           Chemical c2 = Chemical.parse(mol2);
           c2.aromatize();
           try{
        	   c2.generateCoordinates();
           }catch(Exception e){
        	   e.printStackTrace();
           }
           System.out.println(c2.toMol());
           System.out.println(c2.toSmiles());
           Optional<int[]> hit = MolSearcherFactory.create(c).search(c2);
           assertEquals("true",""+hit.isPresent());
       }
       
       
       @Test
       public void ensureAnyBondAndAromatizationInSimpleExampleFromSmartsWorksForSubstructureSearch() throws IOException {
           String mol= "[#7,#8]~C1=CC=CC=C1";
                    
           
	       Chemical c = Chemical.parse(mol);
	       try{
	    	   c.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
	       c.aromatize();
	       System.out.println(c.toMol());
	       System.out.println(c.toSmarts());
//		   IAtomContainer cont = (IAtomContainer) c.getImpl().getWrappedObject();
//		   for(IBond ib : cont.bonds()){
//			   QueryBond qb = ((QueryBond)BondRef.deref(ib));
//			   CdkUtil.navNodes(((QueryBond)qb).getExpression(), (e,l)->{
//					System.out.println(CdkUtil.rep(" ",e) + l.type() + ":" + l.value());
//				},0);
//			   System.out.println(ib.getIndex() + ":" + ib.getOrder() + ":" + ib.isAromatic() + ":" + ib.getClass());
//			   
//		   }
//		   for(IAtom ia : cont.atoms()){
//			   IAtom dr=AtomRef.deref(ia);
//			   if(dr instanceof QueryAtom){
//		    	   QueryAtom qa = (QueryAtom)dr;
//		    	   CdkUtil.navNodes(((QueryAtom)qa).getExpression(), (e,l)->{
//						System.out.println(CdkUtil.rep(" ",e) + l.type() + ":" + l.value());
//					},0);
//			   }
//			   System.out.println(ia.getAtomicNumber() + ":" + ia.getSymbol() + ":" + ia.isAromatic() + ":" + ia.getClass());
//			   
//		   }
	       
	       String mol2="OC1=CC=CC=C1";
	       Chemical c2 = Chemical.parse(mol2);
	       c2.aromatize();
	       try{
	    	   c2.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
	       System.out.println(c2.toMol());
	       System.out.println(c2.toSmiles());
	       Optional<int[]> hit = MolSearcherFactory.create(c).search(c2);
	       assertEquals("true",""+hit.isPresent());
       }
       

       //
      	//String structure="[#7,#8]~C1=c2c3c(OC([#6])(O)C3=O)cc(O)c2=C(O)\\C=C/1";
      	//structureIndexer.add("1", "COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12");
      	//structureIndexer.add("2", "CC1=C2OC(C)(O)C(=O)C2=C3C4=C(C=C(O)C3=C1O)N5C=CC=CC5=N4");
       @Test
       public void ensureAnyBondAndAromatizationInComplexExampleFromSmartsWorksForSubstructureSearch() throws IOException {
           String mol= "[#7,#8]~C1=c2c3c(OC([#6])(O)C3=O)cc(O)c2=C(O)\\C=C/1";
                    
           
	       Chemical c = Chemical.parse(mol);
	       try{
	    	   c.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
	       c.aromatize();
	       System.out.println(c.toMol());
	       System.out.println(c.toSmarts());

	       
	       String mol2="COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12";
	       Chemical c2 = Chemical.parse(mol2);
	       c2.aromatize();
	       try{
	    	   c2.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
	       System.out.println(c2.toMol());
	       System.out.println(c2.toSmiles());
	       Optional<int[]> hit = MolSearcherFactory.create(c).search(c2);
	       assertTrue(hit.isPresent());
       }
       @Test
       public void ensureAnyBondAndAromatizationInComplexExampleFromSmartsWorksForFingerprint() throws IOException {
           String mol= "[#7,#8]~C1=c2c3c(OC([#6])(O)C3=O)cc(O)c2=C(O)\\C=C/1";
                    
           
	       Chemical c = Chemical.parse(mol);
	       try{
	    	   c.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
	       c.aromatize();
	       System.out.println(c.toMol());
	       System.out.println(c.toSmarts());

	       
	       String mol2="COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12";
	       Chemical c2 = Chemical.parse(mol2);
	       c2.aromatize();
	       try{
	    	   c2.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
	       System.out.println(c2.toMol());
	       System.out.println(c2.toSmiles());


	        Fingerprinter fingerPrinterSub =  Fingerprinters.getFingerprinter(FingerprintSpecification.PATH_BASED.create().setLength(512));
	        
	        Fingerprint fp=fingerPrinterSub.computeFingerprint(c);
	        
	        Fingerprint fp2=fingerPrinterSub.computeFingerprint(c2);
	        
	        BitSet bsTemp = BitSet.valueOf(fp.toBitSet().toLongArray());
	        bsTemp.and(fp2.toBitSet());

	        assertEquals(fp.populationCount(), bsTemp.cardinality());
	        
       }
       @Test
	   public void removeAtomThenReAdd() throws Exception{
			Chemical c=Chemical.createFromSmiles("CCCCC");
			Atom remAt=c.getAtom(3);
			List<Bond> removedBonds = new ArrayList<>();
			for(Bond b : remAt.getBonds()){
				removedBonds.add(c.removeBond(b));
			}

			c.removeAtom(remAt);

			c.addAtom(remAt);
			for(Bond b : removedBonds){
				c.addBond(b);
			}
			assertEquals("CCCCC", c.toSmiles());
		}
}
