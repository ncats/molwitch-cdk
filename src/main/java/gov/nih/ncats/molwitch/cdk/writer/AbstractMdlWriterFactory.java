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
import java.util.function.Function;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.IChemObjectWriter;
import org.openscience.cdk.io.listener.PropertiesListener;

import gov.nih.ncats.witch.io.ChemFormat.AromaticAwareChemFormatWriterSpecification;
import gov.nih.ncats.witch.io.ChemFormat.ChemFormatWriterSpecification;
import gov.nih.ncats.witch.io.ChemFormat.KekulizationEncoding;
import gov.nih.ncats.witch.io.ChemFormat.MolFormatSpecification;
import gov.nih.ncats.witch.io.ChemFormat.MolFormatSpecification.Version;
import gov.nih.ncats.witch.spi.ChemicalWriterImpl;
import gov.nih.ncats.witch.spi.ChemicalWriterImplFactory;

public abstract class AbstractMdlWriterFactory implements ChemicalWriterImplFactory{


	protected abstract boolean supportsVersion(Version version);

	@Override
	public ChemicalWriterImpl newInstance(OutputStream out, ChemFormatWriterSpecification spec) throws IOException {
		
		
		 IChemObjectWriter writer = create(out);
		 Function<IAtomContainer, IAtomContainer> adapter = Function.identity();
		 
		 if(spec instanceof AromaticAwareChemFormatWriterSpecification) {
			 Properties customSettings = new Properties();
			 customSettings.setProperty(
			  "WriteAromaticBondTypes", Boolean.toString(((AromaticAwareChemFormatWriterSpecification)spec).getKekulization() == KekulizationEncoding.FORCE_AROMATIC)
			 );
			 PropertiesListener listener =  new PropertiesListener(customSettings);
			 //OptForceWriteAs2DCoordinates - default to false
			 if(spec instanceof MolFormatSpecification) {
				 MolFormatSpecification molSpec = (MolFormatSpecification) spec;
				 
				 switch(molSpec.getCoordinateOptions()) {
				 case FORCE_2D:
					 customSettings.setProperty(
							  "OptForceWriteAs2DCoordinates", Boolean.toString(true)
							 );
					 adapter = c->{
						 //force all coords as 2d
						 boolean mustChange = false;
						 for(IAtom a :c.atoms()) {
							 if(a.getPoint2d() ==null && a.getPoint3d() !=null) {
								 mustChange = true;
							 }
						 }
						 if(mustChange) {
							IAtomContainer clone;
							try {
								clone = c.clone();
							} catch (CloneNotSupportedException e) {
								return c;
							}
							for(IAtom a :clone.atoms()) {
								Point3d points = a.getPoint3d();
								if(points !=null) {
									a.setPoint2d(new Point2d(points.x, points.y));
								}
							}
							return clone;
						 }else {
							 return c;
						 }
					 };
					 break;
				 default:
					 customSettings.setProperty(
							  "OptForceWriteAs2DCoordinates", Boolean.toString(false)
							 );
					 adapter = c->{
						 //force all coords as 2d
						 boolean mustChange = false;
						 for(IAtom a :c.atoms()) {
							 if(a.getPoint3d() ==null && a.getPoint2d() !=null) {
								 mustChange = true;
							 }
						 }
						 if(mustChange) {
							IAtomContainer clone;
							try {
								clone = c.clone();
							} catch (CloneNotSupportedException e) {
								return c;
							}
							for(IAtom a :clone.atoms()) {
								Point2d points = a.getPoint2d();
								if(points !=null) {
									a.setPoint3d(new Point3d(points.x, points.y, 0));
								}
							}
							return clone;
						 }else {
							 return c;
						 }
					 };
					 
					 break;
				 }
			 }
			 writer.addChemObjectIOListener(listener);
		 }
		
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
