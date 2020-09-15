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

import gov.nih.ncats.molwitch.Chemical;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.IOException;

public class CdkUtil {

	public static IAtomContainer toAtomContainer(Chemical chemical){
		return (IAtomContainer) chemical.getImpl().getWrappedObject();
	}
	public static IChemObjectBuilder getChemObjectBuilder() {
		return SilentChemObjectBuilder.getInstance();
	}

	public static IAtomContainer parseSmarts(String smarts) throws CDKException, IOException {
		IAtomContainer container = CdkUtil.getChemObjectBuilder().newAtomContainer();
		Smarts.parse(container, smarts);

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
		QueryAtomPerceptor.percieve(container);
		return container;
	}
}
