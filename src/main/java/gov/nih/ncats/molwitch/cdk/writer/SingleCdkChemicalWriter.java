/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2024.
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
import java.util.List;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.IChemObjectWriter;

import gov.nih.ncats.molwitch.spi.ChemicalImpl;
import org.openscience.cdk.sgroup.Sgroup;
import org.openscience.cdk.sgroup.SgroupKey;

public class SingleCdkChemicalWriter extends CdkChemicalWriter{

	private boolean writtenAlready = false;
	
	public SingleCdkChemicalWriter(IChemObjectWriter writer) {
		super(writer);
	}

	@Override
	public void write(ChemicalImpl impl) throws IOException {
		if(writtenAlready){
			throw new OnlyOneChemicalAllowedException("already wrote a Chemical to this writer");
		}
		super.write(impl);

		writtenAlready=true;
	}

	private static class OnlyOneChemicalAllowedException extends IOException{

		private static final long serialVersionUID = 1L;

		public OnlyOneChemicalAllowedException(String message) {
			super(message);
		}
		
	}
	
}
