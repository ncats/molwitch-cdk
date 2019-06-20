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

import java.io.IOException;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.IChemObjectWriter;

import gov.nih.ncats.witch.spi.ChemicalImpl;
import gov.nih.ncats.witch.spi.ChemicalWriterImpl;
import gov.nih.ncats.molwitch.cdk.CdkChemicalImpl;

public class CdkChemicalWriter implements ChemicalWriterImpl{

	private final IChemObjectWriter writer;

	public CdkChemicalWriter(IChemObjectWriter writer) {
		this.writer = writer;
	}

	@Override
	public void close() throws IOException {
		writer.close();
	}

	@Override
	public void write(ChemicalImpl impl) throws IOException {
		CdkChemicalImpl chem =(CdkChemicalImpl)impl;
		IAtomContainer mol =chem.getContainer();
		try {
			writer.write(mol);
		}catch(Throwable e) {
			throw new IOException("error writing container " + mol.getID(), e);
		}
		
	}

	
	
}
