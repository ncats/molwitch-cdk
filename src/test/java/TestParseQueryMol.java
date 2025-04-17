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

import java.io.IOException;
import java.util.*;

import org.junit.Test;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Bond.BondType;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.fingerprint.Fingerprint;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.FingerprintSpecification;
import gov.nih.ncats.molwitch.inchi.InChiResult;
import gov.nih.ncats.molwitch.search.MolSearcherFactory;

import static org.junit.Assert.*;

public class TestParseQueryMol {

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
        c.aromatize();
//        System.out.println(c.toMol());

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
        System.out.println(mm);
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
        
        assertFalse(smarts.isEmpty());
        
        
        
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
        assertFalse(c2.toSmarts().isEmpty());

        
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

        
        Chemical c2 = Chemical.parseMol(mm);
        
        Optional<int[]> hit = MolSearcherFactory.create(c).get().search(c2);
		

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

           assertFalse(c.toMol().isEmpty());
           assertFalse(c.toSmarts().isEmpty());
           
           
           String mol2="NC1=CC(O)=C(O)C=C1";
           Chemical c2 = Chemical.parse(mol2);
           c2.aromatize();
           try{
        	   c2.generateCoordinates();
           }catch(Exception e){
        	   e.printStackTrace();
           }
		   assertFalse(c2.toMol().isEmpty());
		   assertFalse(c2.toSmiles().isEmpty());
           Optional<int[]> hit = MolSearcherFactory.create(c).get().search(c2);
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
		   assertFalse(c.toMol().isEmpty());
		   assertFalse(c.toSmarts().isEmpty());
           
           
           String mol2="COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12";
           Chemical c2 = Chemical.parse(mol2);
           c2.aromatize();
           try{
        	   c2.generateCoordinates();
           }catch(Exception e){
        	   e.printStackTrace();
           }
		   assertFalse(c2.toMol().isEmpty());
		   assertFalse(c2.toSmiles().isEmpty());
           Optional<int[]> hit = MolSearcherFactory.create(c).get().search(c2);
           assertTrue(hit.isPresent());
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
		   assertFalse(c.toMol().isEmpty());
		   assertFalse(c.toSmarts().isEmpty());

	       
	       String mol2="OC1=CC=CC=C1";
	       Chemical c2 = Chemical.parse(mol2);
	       c2.aromatize();
	       try{
	    	   c2.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
		   assertFalse(c2.toMol().isEmpty());
		   assertFalse(c2.toSmiles().isEmpty());
	       Optional<int[]> hit = MolSearcherFactory.create(c).get().search(c2);
	       assertTrue(hit.isPresent());
       }

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
		   assertFalse(c.toMol().isEmpty());
		   assertFalse(c.toSmarts().isEmpty());

	       
	       String mol2="COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12";
	       Chemical c2 = Chemical.parse(mol2);
	       c2.aromatize();
	       try{
	    	   c2.generateCoordinates();
	       }catch(Exception e){
	    	   e.printStackTrace();
	       }
		   assertFalse(c2.toMol().isEmpty());
		   assertFalse(c2.toSmiles().isEmpty());
	       Optional<int[]> hit = MolSearcherFactory.create(c).get().search(c2);
	       assertTrue(hit.isPresent());
       }
       @Test
       public void ensureAnyBondAndAromatizationInComplexExampleFromSmartsWorksForFingerprint() throws IOException {
           String mol = "[#7,#8]~C1=c2c3c(OC([#6])(O)C3=O)cc(O)c2=C(O)\\C=C/1";

           Chemical c = Chemical.parse(mol);
           try {
               c.generateCoordinates();
           } catch (Exception e) {
               e.printStackTrace();
           }
           c.aromatize();
           assertFalse(c.toMol().isEmpty());
           assertFalse(c.toSmarts().isEmpty());

           String mol2 = "COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12";
           Chemical c2 = Chemical.parse(mol2);
           c2.aromatize();
           try {
               c2.generateCoordinates();
           } catch (Exception e) {
               e.printStackTrace();
           }
           assertFalse(c2.toMol().isEmpty());
           assertFalse(c2.toSmiles().isEmpty());

           Fingerprinter fingerPrinterSub = Fingerprinters.getFingerprinter(
                   FingerprintSpecification.PATH_BASED.create().setLength(512));

           Fingerprint fp = fingerPrinterSub.computeFingerprint(c);

           Fingerprint fp2 = fingerPrinterSub.computeFingerprint(c2);

           double sim = fp.tanimotoSimilarity(fp2);
           OptionalDouble shortSim = fp.tanimotoSimilarityShortCircuit(fp2);
           BitSet bsTemp = fp.toBitSet();
           bsTemp.and(fp2.toBitSet());

           assertEquals(fp.populationCount(), bsTemp.cardinality());
           
           //Depending on how this is interpretted, it might actually be pretty low
           //similarity, and that's okay
           
//           System.out.println(sim);
//           assertTrue(sim > .8D);

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

	   	@Test
	   	public void testQueryAtomMolfileHasWriteAtoms() throws Exception {
	   		Chemical c=Chemical.parse("S(=O)(=O)(O)OC[#6]");

	   		Chemical c2= Chemical.parse(c.toMol());
	   		boolean hasSulfur = c2.atoms().filter(ca->"S".equals(ca.getSymbol())).count()>0;
	       	assertTrue("Simple SMARTS keeps its atom types", hasSulfur);
	   		

	
	   	}
	   	
	  	@Test
	   	public void testQueryAtomMolfileHasWriteAtoms2() throws Exception {
	   		Chemical c=Chemical.parse("S(=O)(=O)(O)OC[#6,#7]");

	   		Chemical c2= Chemical.parse(c.toMol());
	   		boolean hasSulfur = c2.atoms().filter(ca->"S".equals(ca.getSymbol())).count()>0;
	       	assertTrue("Simple SMARTS keeps its atom types", hasSulfur);
			System.out.printf("output molfile: %s\n", c.toMol());
	       	assertTrue("Simple SMARTS keeps atom list", c.toMol().contains("0.0000 L   0"));


	
	   	}

	   	@Test
		public void legacyAtomList() throws Exception{
    	String mol = "\n" +
				"  ACCLDraw10012012192D\n" +
				"\n" +
				" 10 10  3  0  0  0  0  0  0  0999 V2000\n" +
				"   12.8286   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.0309   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5986   -6.6106    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5909   -8.0335    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.8363   -8.7514    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0379   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0714   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.4779   -7.3439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
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
				"M  END\n";

    	Chemical c = Chemical.parse(mol);
    	assertTrue(c.getAtom(5).isQueryAtom());
		}

	@Test
	public void legacyAtomListGetStereoCentersDoentErrorOut() throws Exception{
		String mol = "\n" +
				"  ACCLDraw10012012192D\n" +
				"\n" +
				" 10 10  3  0  0  0  0  0  0  0999 V2000\n" +
				"   12.8286   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.0309   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5986   -6.6106    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5909   -8.0335    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.8363   -8.7514    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0379   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0714   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.4779   -7.3439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
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
				"M  END\n";

		Chemical c = Chemical.parse(mol);
		assertEquals(0, c.getAllStereocenters().size());
	}

	@Test
	public void modernAtomListGetStereoCentersDoentErrorOut() throws Exception{
		String mol = "\n" +
				"  ACCLDraw10012012192D\n" +
				"\n" +
				" 10 10  3  0  0  0  0  0  0  0999 V2000\n" +
				"   12.8286   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.0309   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5986   -6.6106    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5909   -8.0335    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.8363   -8.7514    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0379   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0714   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.4779   -7.3439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
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
				"M  ALS   4  2 F O   N   \n" +
				"M  ALS   5  2 F N   O   \n" +
				"M  ALS   6  2 F N   O   \n" +
				"M  END\n";

		Chemical c = Chemical.parse(mol);
		assertEquals(0, c.getAllStereocenters().size());
	}

	@Test
	public void mixOfModernAtomListGetStereoCentersDoentErrorOut() throws Exception{
		String mol = "\n" +
				"  ACCLDraw10012012192D\n" +
				"\n" +
				" 10 10  3  0  0  0  0  0  0  0999 V2000\n" +
				"   12.8286   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.0309   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5986   -6.6106    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5909   -8.0335    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.8363   -8.7514    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0379   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0714   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.4779   -7.3439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
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
				"M  END\n";

		Chemical c = Chemical.parse(mol);
		assertEquals(0, c.getAllStereocenters().size());
	}

	@Test
	public void clonedMixOfModernAtomListGetStereoCentersDoentErrorOut() throws Exception{
		String mol = "\n" +
				"  ACCLDraw10012012192D\n" +
				"\n" +
				" 10 10  3  0  0  0  0  0  0  0999 V2000\n" +
				"   12.8286   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.0309   -7.3697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5986   -6.6106    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   11.5909   -8.0335    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   12.8363   -8.7514    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   13.2300   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0379   -8.0490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.0714   -6.6260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
				"   14.4779   -7.3439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
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
				"M  END\n";

		Chemical c = Chemical.parse(mol).copy();
		assertEquals(0, c.getAllStereocenters().size());
		assertNotNull(c.getMass());
	}
	   	
}
