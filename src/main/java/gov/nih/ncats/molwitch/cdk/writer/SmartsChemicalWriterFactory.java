/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2023.
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

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smarts.Smarts;

import gov.nih.ncats.molwitch.cdk.CdkUtil;
import gov.nih.ncats.molwitch.io.ChemFormat;
import gov.nih.ncats.molwitch.spi.ChemicalImpl;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImpl;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImplFactory;

/**
 * Created by katzelda on 7/15/19.
 */
public class SmartsChemicalWriterFactory implements ChemicalWriterImplFactory{
    @Override
    public ChemicalWriterImpl newInstance(OutputStream out, ChemFormat.ChemFormatWriterSpecification spec) throws IOException {
        ChemFormat.SmartsFormatSpecification smartSpec = (ChemFormat.SmartsFormatSpecification)spec;

        return new SmartsWriterImpl(out);
    }

    @Override
    public boolean supports(ChemFormat.ChemFormatWriterSpecification spec) {
        return spec instanceof ChemFormat.SmartsFormatSpecification;
    }

    private static class SmartsWriterImpl implements ChemicalWriterImpl{

        private final PrintWriter out;

        public SmartsWriterImpl(OutputStream out) {
            this.out = new PrintWriter(out);
        }
 
        @Override
        public void write(ChemicalImpl chemicalImpl) throws IOException {
        	IAtomContainer ia = CdkUtil.asQueryAtomContainer((IAtomContainer) chemicalImpl.getWrappedObject());
        	
        	out.print(Smarts.generate(ia));
        }

        @Override
        public void close() throws IOException {
            out.close();
        }
    }
}
