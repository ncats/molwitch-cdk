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

import java.io.BufferedReader;
import java.io.IOException;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import gov.nih.ncats.molwitch.cdk.CdkChemical2FactoryImpl.SavedBufferedReader;
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
	GuessResult guessReaderFor(BufferedReader reader) throws IOException{
		reader.mark(Short.MAX_VALUE); //should be big enough...
		
		//a mol file has a header on line 4 and mol2000 and 3000 will say which version
		//any calls to readLine() after we finished reading should return null
		//so multiple calls should be safe
		reader.readLine();
		reader.readLine();
		reader.readLine();
		String header = reader.readLine();
		
		//either way reset so we can parse from the beginning
		reader.reset();
		if(header !=null && (header.indexOf("V2000") >=0 || header.indexOf("V3000") >=0)){
				SavedBufferedReader savedReader = new SavedBufferedReader(new ProgramClearingMol2000Wrapper(reader));
				return new GuessResult(new IdAwareSdfReader(savedReader, SilentChemObjectBuilder.getInstance()), savedReader);
			
		}
		//assume smiles file?
		SavedBufferedReader savedReader = new SavedBufferedReader(reader);
		return new GuessResult(new IteratingSMILESReader( savedReader, SilentChemObjectBuilder.getInstance()), savedReader);
		
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
	
	
	
}
