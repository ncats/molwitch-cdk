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

package gov.nih.ncats.molwitch.cdk.fingerprinters;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import gov.nih.ncats.common.sneak.Sneak;
import gov.nih.ncats.common.util.CachedSupplier;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.ChemicalSource;
import gov.nih.ncats.molwitch.cdk.CdkChemicalImpl;
import gov.nih.ncats.molwitch.cdk.CdkUtil;
import gov.nih.ncats.molwitch.fingerprint.Fingerprint;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;
import gov.nih.ncats.molwitch.spi.FingerprinterImpl;

/**
 * Adapter to wrap a {@link FingerprinterImpl} from the SPI
 * and adapt it to use the {@link Fingerprinter} interface
 * from the API.
 * 
 * @author katzelda
 *
 */
class FingerprinterAdapter implements Fingerprinter {
	private final IFingerprinter delegate;
	private boolean removeQueryAtomsAndBonds=false;
	private boolean forceExplicitH =false;

	public FingerprinterAdapter(IFingerprinter impl) {
		this.delegate = impl;
	}

	public boolean isRemoveQueryAtomsAndBonds() {
		return removeQueryAtomsAndBonds;
	}

	public void setRemoveQueryAtomsAndBonds(boolean removeQueryAtomsAndBonds) {
		this.removeQueryAtomsAndBonds = removeQueryAtomsAndBonds;
	}

	public boolean isForceExplicitH() {
		return forceExplicitH;
	}

	public void setForceExplicitH(boolean forceExplicitH) {
		this.forceExplicitH = forceExplicitH;
	}
	
	private void printContainer(IAtomContainer cln) {
	    
        CdkChemicalImpl cci = new CdkChemicalImpl(cln, (ChemicalSource)null);
        Chemical cq = new Chemical(cci);
        Chemical qq=cq.copy();
        try {
            System.out.println(qq.toMol());
        }catch(Exception e) {
            
        }
	}

	@Override
	public Fingerprint computeFingerprint(Chemical chemical) {
	   
		//cdk doesn't seem to support fingerprints with query atoms/bonds
		//so lets remove them
		IAtomContainer orig = (IAtomContainer)chemical.getImpl().getWrappedObject();

		CachedSupplier<IAtomContainer> cloneCache = CachedSupplier.of(()->{
			try {
			    Chemical cc=chemical.copy();
			    IAtomContainer cln = (IAtomContainer)cc.getImpl().getWrappedObject();
			    try {
		            CdkUtil.removeQueryAtoms(cln);
		            CdkUtil.removeQueryBonds(cln);
		            cln=CdkUtil.getSimplifiedContainer(cln);
			        // The logic here is a little suspect, and may not be necessary,
			        // but the basic idea is to flag all bonds that will
			        // be modified by the aromatic detection and mark them for removal
			        // as they will interfere with fingerprints
			        Map<IBond, String> oldOrder = new HashMap<>();
			        Iterator<IBond> bondIter2 = cln.bonds().iterator();
                    while(bondIter2.hasNext()){
                        IBond next = bondIter2.next();
                        oldOrder.put(next, next.getOrder() + "" + next.isAromatic());
                    }
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(cln);
		            Aromaticity.cdkLegacy().apply(cln);
		            CdkUtil.removeQueryBonds(cln);
		           
		            Iterator<IBond> bondIter = cln.bonds().iterator();
		            while(bondIter.hasNext()){
		                IBond next = bondIter.next();
		                if(next.getOrder() ==null || next.getOrder() == Order.UNSET ||
		                        !(next.getOrder() + "" + next.isAromatic()).equals(oldOrder.get(next))){
		                    next.setOrder(Order.UNSET);
		                }
		            }
		        }catch(Exception e) {
		            //swallow issue
//		            e.printStackTrace();		            
		        }
				return cln;
			}catch(Exception t){
				return Sneak.sneakyThrow(t);
			}
		});
		IAtomContainer containerToFingerprint = orig;

		if(removeQueryAtomsAndBonds){
			containerToFingerprint = cloneCache.get();
			CdkUtil.removeQueryAtoms(containerToFingerprint);
			CdkUtil.removeQueryBonds(containerToFingerprint);
		}
		//This whole idea here is a bit suspicious
		if(forceExplicitH){
			containerToFingerprint = cloneCache.get();
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(containerToFingerprint);
		}
		
		IAtomContainer iac= CdkUtil.getUsableFormOfAtomContainer(containerToFingerprint);
		
		// TODO:
		// IF IAtomContainer is a query container, it simply doesn't produce
		// ANY fingerprint, and something needs to be done.

		if(iac instanceof QueryAtomContainer) {
		    iac = CdkUtil.getSimplifiedContainer(iac);
		}
		
		try {
		    IBitFingerprint bitFingerprint = delegate.getBitFingerprint(iac);
		    BitSet bs = bitFingerprint.asBitSet();
		    return new Fingerprint(bs, (int) bitFingerprint.size());
		} catch (CDKException e) {
		    throw new RuntimeException(e);
		}
	}
}