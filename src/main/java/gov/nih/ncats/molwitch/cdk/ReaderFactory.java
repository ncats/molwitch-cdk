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

import java.io.*;
import java.util.NoSuchElementException;
import java.util.stream.Stream;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.formats.*;
import org.openscience.cdk.io.iterator.DefaultIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import gov.nih.ncats.molwitch.cdk.CdkChemical2FactoryImpl.SavedBufferedReader;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * A Factory class for making {@link IIteratingChemObjectReader}s.
 * For some reason, CDK (as of 1.5.13) doesn't have a factory class
 * that can guess the file format and return the appropriate
 * {@link IIteratingChemObjectReader} for it, but it does have one
 * for the non-iterating, read everything into memory Reader.
 *
 * <p>
 * This is a very basic factory that can only determine
 * a file encoded in SMILES format or MOL/SDF format.
 * If the file is compressed, it must be uncompressed before it is
 * given to this factory.
 *
 * @author katzelda
 *
 */
public class ReaderFactory {


	/**
	 * Get the file format for the contents of the given {@link BufferedReader}
	 * and return the appropriate {@link IIteratingChemObjectReader} for that format.
	 *
	 * @param reader the {@link BufferedReader} to parse; can not be null
	 * and must support mark/reset.
	 *
	 * @return a new {@link IIteratingChemObjectReader} that can parse
	 * the input.
	 * @throws IOException if there is a problem reading the input data.
	 */
	public static GuessResult guessReaderFor(BufferedReader reader) throws IOException{
		//CDK reader gets confused by mol files with a name in the first line because
		//depending on the format of the name it can guess the wrong format
		//so check for mol file first and if it's not then use cdk to guess
		IResourceFormat format=null;
		reader.mark(2000);
		reader.readLine();
		reader.readLine();
		reader.readLine();
		String molLine = reader.readLine();
		if(molLine !=null) {
			if (molLine.endsWith("V2000")) {
				format = MDLV2000Format.getInstance();
			} else if (molLine.endsWith("V3000")) {
				format = MDLV3000Format.getInstance();
			}
		}
		if(format ==null) {
			try {
				reader.reset();
				format = new FormatFactory().guessFormat(reader);
			} catch (Throwable e) {
				e.printStackTrace();
				//TODO should probably make a smiles or smarts reader...

				format = SMILESFormat.getInstance();
			}
		}
		if(format ==null){
			//default to smiles?
			format = SMILESFormat.getInstance();
		}

		reader.reset();
		return create(reader, format);

	}

	public static GuessResult create(BufferedReader reader, IResourceFormat format) throws IOException {
		if(format instanceof MDLV2000Format || format instanceof MDLV3000Format){
			SavedBufferedReader savedReader = new SavedBufferedReader(new ProgramClearingMol2000Wrapper(reader));
			return new GuessResult(new IdAwareSdfReader(savedReader, SilentChemObjectBuilder.getInstance()), savedReader);

		}

		if(format instanceof SMILESFormat){
			SavedBufferedReader savedReader = new SavedBufferedReader(reader);
			return new GuessResult(new IteratingSMILESReader( savedReader, SilentChemObjectBuilder.getInstance()), savedReader);

		}
		if(format instanceof SMARTSFormat) {
			SavedBufferedReader savedReader = new SavedBufferedReader(reader);
			return new GuessResult(new IteratingSmartsReader(savedReader, SilentChemObjectBuilder.getInstance()), savedReader);

		}
		if(format instanceof SDFFormat){
			SavedBufferedReader savedReader = new SavedBufferedReader(new ProgramClearingMol2000Wrapper(reader));
			return new GuessResult(new IdAwareSdfReader(savedReader, SilentChemObjectBuilder.getInstance()), savedReader);

		}
		throw new IOException ("not configured to parse " + format.getFormatName());
	}

	static class GuessResult{
		public final IIteratingChemObjectReader<IAtomContainer> cdkReader;
		public final SavedBufferedReader savedBufferedReader;
		public GuessResult(IIteratingChemObjectReader<IAtomContainer> cdkReader,
						   SavedBufferedReader savedBufferedReader) {
			this.cdkReader = cdkReader;
			this.savedBufferedReader = savedBufferedReader;
		}


	}

	static class IteratingSmartsReader extends DefaultIteratingChemObjectReader<IAtomContainer>{

		private BufferedReader input;
		private IAtomContainer nextMolecule;
		private boolean nextAvailableIsKnown, hasNext;
		private final IChemObjectBuilder builder;

		public IteratingSmartsReader(Reader in, IChemObjectBuilder builder) {

			this.setReader(in);
			this.builder = builder;
		}

		public IteratingSmartsReader(InputStream in, IChemObjectBuilder builder) {
			this(new InputStreamReader(in), builder);
		}
		public boolean hasNext() {
			if (!this.nextAvailableIsKnown) {
				this.hasNext = false;

				try {
					String line = this.input.readLine();
					if (line == null) {
						this.nextAvailableIsKnown = true;
						return false;
					}

					this.hasNext = true;
					IAtomContainer container = builder.newAtomContainer();
					Smarts.parse(container, line);

					AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
					QueryAtomPerceptor.percieve(container);
					this.nextMolecule = container;
				} catch (Exception var3) {
					this.hasNext = false;
				}

				if (!this.hasNext) {
					this.nextMolecule = null;
				}

				this.nextAvailableIsKnown = true;
			}

			return this.hasNext;
		}
		public IAtomContainer next() {
			if (!this.nextAvailableIsKnown) {
				this.hasNext();
			}

			this.nextAvailableIsKnown = false;
			if (!this.hasNext) {
				throw new NoSuchElementException();
			} else {
				return this.nextMolecule;
			}
		}

		public void close() throws IOException {
			if (this.input != null) {
				this.input.close();
			}

		}


		public void setReader(Reader reader) {
			if (reader instanceof BufferedReader) {
				this.input = (BufferedReader)reader;
			} else {
				this.input = new BufferedReader(reader);
			}

			this.nextMolecule = null;
			this.nextAvailableIsKnown = false;
			this.hasNext = false;
		}

		public void setReader(InputStream reader) {
			this.setReader((Reader)(new InputStreamReader(reader)));
		}

		@Override
		public IResourceFormat getFormat() {
			return SMARTSFormat.getInstance();
		}

	}



}
