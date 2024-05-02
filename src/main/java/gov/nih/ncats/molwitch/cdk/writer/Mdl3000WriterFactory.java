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
import java.io.OutputStream;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.IChemObjectWriter;
import org.openscience.cdk.io.MDLV3000Writer;

import gov.nih.ncats.molwitch.io.ChemFormat.MolFormatSpecification.Version;

public class Mdl3000WriterFactory extends AbstractMdlWriterFactory{

	


	@Override
	protected boolean supportsVersion(Version version) {
		return version == Version.V3000;
	}

	@Override
	protected IChemObjectWriter create(OutputStream out) throws IOException {
		try {
			return new MDLV3000Writer(out);
		} catch (CDKException e) {
			throw new IOException("error creating mol V3000 writer", e);
		}
	}

}
