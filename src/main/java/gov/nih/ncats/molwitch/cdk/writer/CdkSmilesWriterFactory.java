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

package gov.nih.ncats.molwitch.cdk.writer;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Properties;
import java.util.function.Function;

import gov.nih.ncats.molwitch.cdk.CdkUtil;
import org.jooq.lambda.Unchecked;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import gov.nih.ncats.molwitch.io.ChemFormat.ChemFormatWriterSpecification;
import gov.nih.ncats.molwitch.io.ChemFormat.HydrogenEncoding;
import gov.nih.ncats.molwitch.io.ChemFormat.KekulizationEncoding;
import gov.nih.ncats.molwitch.io.ChemFormat.SmilesFormatWriterSpecification;
import gov.nih.ncats.molwitch.spi.ChemicalImpl;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImpl;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImplFactory;

public class CdkSmilesWriterFactory implements ChemicalWriterImplFactory{



	@Override
	public boolean supports(ChemFormatWriterSpecification spec) {
		return spec instanceof SmilesFormatWriterSpecification;
	}

	@Override
	public ChemicalWriterImpl newInstance(OutputStream out, ChemFormatWriterSpecification spec) throws IOException {
		SmilesFormatWriterSpecification smilesSpec = (SmilesFormatWriterSpecification) spec;
		int options =SmiFlavor.Generic;
		if(smilesSpec.getCanonization() == SmilesFormatWriterSpecification.CanonicalizationEncoding.CANONICAL) {
			options = SmiFlavor.Canonical;
		}
		if(smilesSpec.getEncodeStereo() == SmilesFormatWriterSpecification.StereoEncoding.INCLUDE_STEREO) {
			options |= SmiFlavor.Stereo;
		}
		if(smilesSpec.getKekulization() == KekulizationEncoding.FORCE_AROMATIC) {
			options |= SmiFlavor.UseAromaticSymbols;
		}
		//include isotope information
		options |= SmiFlavor.AtomicMass;
		return new CdkSmilesWriter(out, options, smilesSpec.getHydrogenEncoding(), smilesSpec.getKekulization());
	}

	
	private static class CdkSmilesWriter implements ChemicalWriterImpl{

		private final PrintWriter out;
		private final SmilesGenerator sg;
		private final Function<IAtomContainer, IAtomContainer> modificationFunction;
		
		private static boolean hasImplicitH(IAtomContainer container) {
			for(IAtom atom: container.atoms()) {
				 Integer implicitNum = atom.getImplicitHydrogenCount();
		            if(implicitNum !=null && implicitNum.intValue() >0){
		            	return true;
		            }
			}
			return false;
		}
		public CdkSmilesWriter(OutputStream out, int flavor, 
				HydrogenEncoding hydrogenEncoding, KekulizationEncoding aromaticEncoding) {
			this.out = new PrintWriter(out);
			sg = new SmilesGenerator(flavor);
			//there is a lot of code duplication here
			//this is mostly because I'd rather have a little code copy + paste
			//and have an easy way to make sure I clone the least number of times
			if (hydrogenEncoding == HydrogenEncoding.MAKE_EXPLICIT) {
				modificationFunction = Unchecked.function(container -> {
					if (hasImplicitH(container)) {
						IAtomContainer copy = container.clone();

						AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);

						if(aromaticEncoding == KekulizationEncoding.FORCE_AROMATIC) {
                            CdkUtil.aromatize(copy);
						}
						return copy;
					}
					return container;

				});
			} else if (hydrogenEncoding == HydrogenEncoding.MAKE_IMPLICIT) {
				modificationFunction = Unchecked.function(container -> {

					IAtomContainer copy = container.clone();
					AtomContainerManipulator.suppressHydrogens(copy);
					if(aromaticEncoding == KekulizationEncoding.FORCE_AROMATIC) {
						boolean alreadyMarkedAromatic=false;
						for(IBond b: copy.bonds()) {
							if(b.isAromatic()) {
								alreadyMarkedAromatic=true;
								break;
							}
						}
						if(!alreadyMarkedAromatic) {
                            CdkUtil.aromatize(copy);
						}
					}
					return copy;

				});
			}else if(aromaticEncoding == KekulizationEncoding.FORCE_AROMATIC) {
				modificationFunction = Unchecked.function(container ->{
				boolean alreadyMarkedAromatic=false;
					for(IBond a: container.bonds()) {
						if(a.isAromatic()) {
							alreadyMarkedAromatic=true;
							break;
						}
					}
					if(alreadyMarkedAromatic) {
						return container;
					}

						IAtomContainer copy = container.clone();

					 CdkUtil.aromatize(copy);
					 return copy;
					
					
				});
			
			}else if(aromaticEncoding == KekulizationEncoding.KEKULE){
				modificationFunction = container ->{
					try{
						return  CdkUtil.kekulizeIfNeeded(container,true);
					}catch(Exception e){
						return container;
					}
				};
			}else {

				modificationFunction = Function.identity();
			}
		}

		@Override
		public void close() throws IOException {
			out.close();
			
		}

		@Override
		public void write(ChemicalImpl chemicalImpl) throws IOException {
			//smiles writer has problems with empty
			if(chemicalImpl.getAtomCount() ==0){
				out.println("");
				return;
			}
			
			//TODO: having unset implicit H count makes this fail
			//but it really probably shouldn't.
			
			try {
				out.println(sg.create(modificationFunction.apply((IAtomContainer)chemicalImpl.getWrappedObject())));
			} catch (CDKException e) {
				throw new IOException("error writing out smiles for " + chemicalImpl.getName(), e);
			}
			
		}
		
	}

}
