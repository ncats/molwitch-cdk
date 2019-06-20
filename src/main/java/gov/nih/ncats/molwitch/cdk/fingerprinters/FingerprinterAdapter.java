/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2019.
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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.openscience.cdk.AtomRef;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryBond;

import gov.nih.ncats.witch.Chemical;
import gov.nih.ncats.witch.fingerprint.Fingerprint;
import gov.nih.ncats.witch.fingerprint.Fingerprinter;
import gov.nih.ncats.witch.spi.FingerprinterImpl;

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
	
	
	public FingerprinterAdapter(IFingerprinter impl) {
		this.delegate = impl;
	}


	@Override
	public Fingerprint computeFingerprint(Chemical chemical) {
		//cdk doesn't seem to support fingerprints with query atoms/bonds
		//so lets remove them
		IAtomContainer orig = (IAtomContainer)chemical.getImpl().getWrappedObject();
//		System.out.println("=====  " + orig.getAtomCount());
		
		
		IAtomContainer container=orig;
		
		try {
			IAtomContainer clone = orig.clone();
			//iterator.remove() doesn't seem to work
			List<Integer> indexesToRemove = new ArrayList<>();
			int i=0;
			for(IAtom next : clone.atoms()){
				
//				System.out.println("atom  = " + next + " sym = " + next.getSymbol());
				if(next.getSymbol() ==null){
					indexesToRemove.add(Integer.valueOf(i));
				}
				i++;
			}
//			System.out.println(indexesToRemove);
			//iterate in reverse to not mess up atom order of deletions
			for(int k= indexesToRemove.size() -1; k >=0; k-- ){
				clone.removeAtom(k);
			}
			Iterator<IBond> bondIter = clone.bonds().iterator();
			while(bondIter.hasNext()){
				
				IBond next = bondIter.next();
				if(next.getOrder() ==null){
//					System.out.println("removing bond " + next);
					bondIter.remove();
				}
			}
//			System.out.println("# atoms left =" + clone.getAtomCount());
			container = clone;
		} catch (CloneNotSupportedException e1) {
			//shouldn't happen
			throw new IllegalStateException(e1);
		}
		
		
		try {
			return new Fingerprint(delegate.getBitFingerprint(container).asBitSet());
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
	}
}