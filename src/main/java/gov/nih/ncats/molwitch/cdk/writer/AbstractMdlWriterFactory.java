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

import gov.nih.ncats.molwitch.io.ChemFormat.ChemFormatWriterSpecification;
import gov.nih.ncats.molwitch.io.ChemFormat.MolFormatSpecification;
import gov.nih.ncats.molwitch.io.ChemFormat.MolFormatSpecification.Version;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImpl;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImplFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.IChemObjectWriter;
import org.openscience.cdk.io.listener.PropertiesListener;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Properties;
import java.util.function.Function;

public abstract class AbstractMdlWriterFactory implements ChemicalWriterImplFactory{


	protected abstract boolean supportsVersion(Version version);

	@Override
	public ChemicalWriterImpl newInstance(OutputStream out, ChemFormatWriterSpecification spec) throws IOException {
		
		
		IChemObjectWriter writer = create(out);
		Properties customSettings = new Properties();
		Function<IAtomContainer, IAtomContainer> adapter = CtabWriterUtil.handleMolSpec( (MolFormatSpecification) spec, customSettings);
		PropertiesListener listener = new PropertiesListener(customSettings);

		writer.addChemObjectIOListener(listener);

		
		return new SingleCdkChemicalWriter(ChemObjectWriterAdapter.create(writer, adapter));
	}


	@Override
	public boolean supports(ChemFormatWriterSpecification spec) {
		if(! MolFormatSpecification.NAME.equalsIgnoreCase(spec.getFormatName())) {
			return false;
		}
		if(spec instanceof MolFormatSpecification) {
			return supportsVersion(((MolFormatSpecification)spec).getVersion());
		}
		return false;
	}
	


	protected abstract IChemObjectWriter create(OutputStream out)  throws IOException ;

}
