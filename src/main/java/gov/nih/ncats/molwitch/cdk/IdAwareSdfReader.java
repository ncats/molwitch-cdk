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

package gov.nih.ncats.molwitch.cdk;

import java.io.BufferedReader;
import java.io.IOException;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
/**
 * And SDF/ Mol parser file parser
 * that uses the first line in the header record
 * as the ID.  The CDK SDF reader doesn't parse the molecule
 * names from the header.  
 * @author katzelda
 *
 */
public class IdAwareSdfReader extends IteratingSDFReader{

	private boolean alreadyCalledHasNext;
	private BufferedReader in;
	private String nextId;
	
	public IdAwareSdfReader(BufferedReader in, IChemObjectBuilder builder) {
		super(in, builder);
		this.in= in;
	}

	@Override
	public boolean hasNext() {
		if(!alreadyCalledHasNext){
			try {
				in.mark(2048);
				nextId = in.readLine();
			} catch (IOException e) {
				//ignore
			}finally{
				try {
					in.reset();
				} catch (IOException e) {
					//ignore
				}
			}
		}
		boolean next = super.hasNext();
		alreadyCalledHasNext = true;
		return next;
	}

	@Override
	public IAtomContainer next() {
		IAtomContainer ret= super.next();
		ret.setID(nextId);
		//some cdk objects like mol writers look for the title property
		ret.setProperty(CDKConstants.TITLE, nextId);
		
		alreadyCalledHasNext = false;
		return ret;
	}
	
	

}
