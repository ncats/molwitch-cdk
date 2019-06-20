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
import java.io.OutputStream;
import java.util.Properties;

import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.listener.PropertiesListener;

import gov.nih.ncats.witch.io.ChemFormat.ChemFormatWriterSpecification;
import gov.nih.ncats.witch.io.ChemFormat.KekulizationEncoding;
import gov.nih.ncats.witch.io.ChemFormat.MolFormatSpecification.Version;
import gov.nih.ncats.witch.io.ChemFormat.SdfFormatSpecification;
import gov.nih.ncats.witch.spi.ChemicalWriterImpl;
import gov.nih.ncats.witch.spi.ChemicalWriterImplFactory;

public class SdfWriterFactory implements ChemicalWriterImplFactory{

	
	@Override
	public ChemicalWriterImpl newInstance(OutputStream out, ChemFormatWriterSpecification spec) throws IOException {
		SdfFormatSpecification sdfSpec = (SdfFormatSpecification) spec;
		SDFWriter writer = new SDFWriter(out);
		 if(sdfSpec.getMolSpec().getVersion() == Version.V3000) {
			 writer.setAlwaysV3000(true);
		 }
		 if(sdfSpec.getMolSpec().getKekulization() == KekulizationEncoding.FORCE_AROMATIC) {
			 Properties customSettings = new Properties();
			 customSettings.setProperty("WriteAromaticBondTypes", Boolean.TRUE.toString());
			 writer.addChemObjectIOListener( new PropertiesListener(customSettings));
			 
		 }
		 return new CdkChemicalWriter(writer);
	}

	@Override
	public boolean supports(ChemFormatWriterSpecification spec) {
		return SdfFormatSpecification.NAME.equals(spec.getFormatName()) && spec instanceof SdfFormatSpecification;
	}



}
