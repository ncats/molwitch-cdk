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

package gov.nih.ncats.molwitch.cdk.writer;

import gov.nih.ncats.molwitch.io.ChemFormat;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import gov.nih.ncats.molwitch.io.ChemFormat.ChemFormatWriterSpecification;
import gov.nih.ncats.molwitch.io.WriterOptions;

public class CdkWriterOptionUtils {

	
	
	public static IAtomContainer createManipulatedContainer(ChemFormatWriterSpecification spec, IAtomContainer container) {
		return createManipulatedContainer(spec, container, false);
	}
	public static IAtomContainer createManipulatedContainer(ChemFormatWriterSpecification spec, IAtomContainer container, boolean clonedContainer) {
		IAtomContainer ret = container;
		if(spec instanceof ChemFormat.HydrogenAwareChemFormatWriterSpecification){
		//order of options in specific order to limit the number of clones
		//DON'T CHANGE ORDER UNLESS YOU KNOW WHAT YOU'RE DOING
		switch( ((ChemFormat.HydrogenAwareChemFormatWriterSpecification)spec).getHydrogenEncoding()) {
			case AS_IS:
				break;
			case MAKE_IMPLICIT:
				ret = AtomContainerManipulator.removeHydrogens(ret);
				break;
			case MAKE_EXPLICIT: {
				if (!clonedContainer) {
					try {
						ret = ret.clone();
					} catch (CloneNotSupportedException e) {
						throw new IllegalStateException("could not make clone", e);
					}
				}
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(container);
				break;
			}
		}
							
		}
		return ret;
	}
}
