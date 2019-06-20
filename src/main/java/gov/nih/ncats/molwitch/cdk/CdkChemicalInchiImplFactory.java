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

import java.io.IOException;
import java.util.Collections;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.ChemicalBuilder;
import gov.nih.ncats.molwitch.ChemicalSource.Type;
import gov.nih.ncats.molwitch.inchi.InChiResult;
import gov.nih.ncats.molwitch.inchi.InChiResult.Status;
import gov.nih.ncats.molwitch.internal.source.StringSource;
import gov.nih.ncats.molwitch.spi.InchiImplFactory;
import net.sf.jniinchi.INCHI_RET;

public class CdkChemicalInchiImplFactory implements InchiImplFactory{

	private final InChIGeneratorFactory factory;
	
	public CdkChemicalInchiImplFactory() {
		try {
			factory = InChIGeneratorFactory.getInstance();
		} catch (CDKException e) {
			throw new IllegalStateException("could not initialize Inchi generator", e);
		}
		
	}
	@Override
	public InChiResult asStdInchi(Chemical chemical, boolean trustCoordinates) throws IOException {
		
		try {
			//need to pass list options (even empty) to get AuxInfo...
			InChIGenerator gen = factory.getInChIGenerator((IAtomContainer)chemical.getImpl().getWrappedObject(), Collections.emptyList());
		
			InChiResult.Status status = toChemkitStatus(gen.getReturnStatus());
			
			return new InChiResult.Builder(status)
						.setAuxInfo(gen.getAuxInfo())
						.setInchi(gen.getInchi())
						.setKey(gen.getInchiKey())
						.setMessage(gen.getLog())
						.build();
			
		
		} catch (CDKException e) {
			throw new IOException("error computing Inchi", e);
		} 
	}

	private InChiResult.Status toChemkitStatus(INCHI_RET returnStatus) {
		if(returnStatus == INCHI_RET.OKAY) {
			return Status.VALID;
		}
		if(returnStatus == INCHI_RET.WARNING) {
			return Status.WARNING;
		}
		return Status.ERROR;
	}
	@Override
	public Chemical parseInchi(String inchi) throws IOException {
		try {
			InChIToStructure toStruc= factory.getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());
		
			System.out.println(toStruc.getReturnStatus() + "  " + toStruc.getMessage());
			if(toStruc.getReturnStatus() == INCHI_RET.OKAY || toStruc.getReturnStatus() == INCHI_RET.WARNING) {
				
				return ChemicalBuilder._fromImpl(new CdkChemicalImpl( toStruc.getAtomContainer(), new StringSource(inchi, Type.INCHI)))
								.build();
							
			}
			throw new IOException("error trying to parse inchi '"+ inchi + "' : " + toStruc.getMessage());
		
		} catch (CDKException e) {
			throw new IOException("error trying to parse inchi '"+ inchi + "'", e);
		}
	}

}
