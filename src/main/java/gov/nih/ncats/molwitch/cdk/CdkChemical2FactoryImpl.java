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
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.function.Consumer;
import java.util.stream.Stream;

import gov.nih.ncats.molwitch.ChemicalSource;
import gov.nih.ncats.molwitch.SmartsSource;
import gov.nih.ncats.molwitch.ChemicalSource.CommonProperties;
import gov.nih.ncats.molwitch.SmilesSource;
import gov.nih.ncats.molwitch.io.ChemFormat.MolFormatSpecification;
import gov.nih.ncats.molwitch.io.ChemFormat.SdfFormatSpecification;
import gov.nih.ncats.molwitch.io.ChemFormat.SmilesFormatWriterSpecification;
import gov.nih.ncats.molwitch.internal.source.MolStringSource;

import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import gov.nih.ncats.molwitch.spi.ChemicalImpl;
import gov.nih.ncats.molwitch.spi.ChemicalImplFactory;
import gov.nih.ncats.molwitch.spi.ChemicalImplReader;
import gov.nih.ncats.common.io.InputStreamSupplier;

public class CdkChemical2FactoryImpl implements ChemicalImplFactory{

	private SmilesParser preserveAromaticSmilesParser;
	private SmilesParser kekuleSmilesParser;
	
	public CdkChemical2FactoryImpl(){
		
		 preserveAromaticSmilesParser = new SmilesParser(CdkUtil.getChemObjectBuilder());
		 //make aromatic
		 preserveAromaticSmilesParser.kekulise(false);
		 kekuleSmilesParser = new SmilesParser(CdkUtil.getChemObjectBuilder());
	}
	
	@Override
	public ChemicalImpl createNewEmptyChemical() {
		return new CdkChemicalImpl(CdkUtil.getChemObjectBuilder().newAtomContainer(), (ChemicalSource)null);
	}

	private IAtomContainer createIAtomContainerFrom(String smiles) throws Exception {
//		String fixedSmiles;
//		if(smiles.contains("[#")){
//			 Isotopes isotopes = Isotopes.getInstance();
//		//Jchem can understand things like [#6] as Carbon
//			 Pattern atomNumPattern = Pattern.compile("\\[#(\\d+)\\]");
//			 Matcher m = atomNumPattern.matcher(smiles);
//			 StringBuilder builder = new StringBuilder();
//			 int prev=-1;
//			 while(m.find()){
//				 builder.append(smiles.substring(prev+1, m.start(1)-2));
//				
//				 builder.append(isotopes.getElementSymbol(Integer.parseInt(m.group(1))));
//				
//				 prev = m.end(1);
//			 }
//			 //append remaining
//			 builder.append(smiles.substring(prev+1));
//			fixedSmiles = builder.toString();
//		}else{
//			fixedSmiles = smiles;
//		}
		IAtomContainer mol = tryCreate(smiles);
//		 addImplicitHydrogens(mol);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);





		StructureDiagramGenerator coordinateGenerator = new StructureDiagramGenerator(mol);
		coordinateGenerator.generateExperimentalCoordinates();
		//coordinateGenerator.generateCoordinates(coordinateGenerator.DEFAULT_BOND_VECTOR, true, true);
		mol = coordinateGenerator.getMolecule();


		return mol;
	}
	
	private IAtomContainer tryCreate(String smiles) throws InvalidSmilesException{
		try {
			return kekuleSmilesParser.parseSmiles(smiles);
		} catch (InvalidSmilesException e) {
//			e.printStackTrace();
			
			return preserveAromaticSmilesParser.parseSmiles(smiles);
		}
		
	}
	
	 /**
     * Copied from CDK CDKTestCase.java adds implicit hydrogens to given container.
	  * @param container the {@link IAtomContainer} to add the hydrogens to.
	  *
	  * @throws Exception if any problems are encountered adding the hydrogens.
     */
    protected void addImplicitHydrogens(IAtomContainer container) throws Exception {
        CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(container.getBuilder());
        int atomCount = container.getAtomCount();
        String[] originalAtomTypeNames = new String[atomCount];
        for (int i = 0; i < atomCount; i++) {
            IAtom atom = container.getAtom(i);
            IAtomType type = matcher.findMatchingAtomType(container, atom);
            originalAtomTypeNames[i] = atom.getAtomTypeName();
            atom.setAtomTypeName(type.getAtomTypeName());
        }
        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(container.getBuilder());
        hAdder.addImplicitHydrogens(container);
        // reset to the original atom types
        for (int i = 0; i < atomCount; i++) {
            IAtom atom = container.getAtom(i);
            atom.setAtomTypeName(originalAtomTypeNames[i]);
        }
    }
	@Override
	public ChemicalImpl createFromSmiles(String smiles) throws IOException {
		
		try {
			IAtomContainer container = createIAtomContainerFrom(smiles);


			CdkChemicalImpl chem=  new CdkChemicalImpl(container, new SmilesSource(smiles));
			
			return chem;
		} catch (Exception e) {
			throw new IOException("error parsing smiles : " + smiles, e);
		}
	}
	
	

	@Override
	public ChemicalImpl createFromSmarts(String smarts) throws IOException{
		try{
			IAtomContainer container = CdkUtil.getChemObjectBuilder().newAtomContainer();
			Smarts.parse(container, smarts);
			
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
			QueryAtomPerceptor.percieve(container);
			return new CdkChemicalImpl(container, new SmartsSource(smarts));
		}catch(Exception e){
			e.printStackTrace();
			throw new IOException("error parsing smarts : '" + smarts + "' reason = " + Smarts.getLastErrorMesg(), e);
		}
	}

	@Override
	public ChemicalImplReader create(byte[] molBytes, int start, int length) throws IOException {
		return create(new ByteArrayInputStream(molBytes, start, length));
	}

	
	@Override
	public ChemicalImplReader create(String format, InputStreamSupplier in) throws IOException {
		return create(in.get());
	}

	@Override
	public ChemicalImplReader create(InputStreamSupplier in) throws IOException {
		return create(in.get());
	}

	@Override
	public ChemicalImplReader create(String format, InputStream in) throws IOException {
		return create(in);
	}

	@Override
	public ChemicalImplReader create(String format, File file) throws IOException {
		return create(file);
	}

	@Override
	public boolean supports(String format) {
		return SmilesFormatWriterSpecification.NAME.equalsIgnoreCase(format)
				|| MolFormatSpecification.NAME.equalsIgnoreCase(format)
				|| SdfFormatSpecification.NAME.equalsIgnoreCase(format)
				;
	}

	@Override
	public ChemicalImplReader create(File file) throws IOException {
		return createFrom(new InputStreamReader(InputStreamSupplier.forFile(file).get()),
				source -> {
					source.getProperties().put(CommonProperties.Filename, file.getName());
					source.getProperties().put(CommonProperties.Filepath, file.getAbsolutePath());
					source.getProperties().put(CommonProperties.Filesize, Long.toString(file.length()));
					});
	}
	
	private ChemicalImplReader createFrom(Reader reader, Consumer<ChemicalSource> sourceConsumer) throws IOException {
		//needs to be buffered because the reader factory
		//tries to go back and re-read!
		
		ReaderFactory readerFactory = new ReaderFactory();
		ReaderFactory.GuessResult guessedReader = readerFactory.guessReaderFor(new BufferedReader(reader));
		return new CdkChemicalImplReader(guessedReader.cdkReader, guessedReader.savedBufferedReader, sourceConsumer);
	}
	
	@Override
	public ChemicalImplReader create(InputStream in) throws IOException {
		return createFrom(new InputStreamReader(in), null);
	}
	
	private static class CdkChemicalImplReader implements ChemicalImplReader{
		private final IIteratingChemObjectReader<IAtomContainer> iter;
		private final SavedBufferedReader savedReader;
		private ChemicalSource.Type type;
		boolean alreadyReadFirstRecord=false;
		private Consumer<ChemicalSource> sourceConsumer;
		public CdkChemicalImplReader(IIteratingChemObjectReader<IAtomContainer> iter, 
				SavedBufferedReader savedReader,
				Consumer<ChemicalSource> sourceConsumer) {
			this.iter = iter;
			this.savedReader = savedReader;
			this.sourceConsumer = sourceConsumer;
			
			this.type = iter instanceof IteratingSDFReader ? ChemicalSource.Type.SDF : ChemicalSource.Type.SMILES;
		}
		@Override
		public void close() throws IOException {
			iter.close();
		}

		@Override
		public ChemicalImpl read() throws IOException {
			
			if(iter.hasNext()){
				try {
				ChemicalImpl impl= new CdkChemicalImpl(iter.next(),()->{
					String data = savedReader.getBufferedLines();
					
					if(type == ChemicalSource.Type.SMILES){
						//only trim smiles
						//mol +sd files white space in header matters
						return new SmilesSource(data.trim());
					}
					
					if(!alreadyReadFirstRecord) {
							
						if(data.contains("$$$$")){
							//I guess it's an SDF?
							type = ChemicalSource.Type.SDF;
						}else{
							type = ChemicalSource.Type.MOL;
						}
					}
					alreadyReadFirstRecord=true;
					return new MolStringSource(data, type);
				});
				savedReader.resetBuffer();
				
				if(sourceConsumer !=null) {
					ChemicalSource source = impl.getSource();
					if(source !=null) {
						sourceConsumer.accept(source);
					}
				}
				return impl;
				}catch(Throwable e) {
					e.printStackTrace();
					System.out.println("problem record\n=======\n" + savedReader.getBufferedLines() +"\n=======");
					throw e;
				}
			}
			
			return null;
		}
		
		
	}

	static class SavedBufferedReader extends BufferedReader{

		private final StringBuilder buffer = new StringBuilder(2048);
		private final String NEW_LINE = System.lineSeparator();
		
		private int resetPosition = -1;
		BufferedReader reader;
		public SavedBufferedReader(BufferedReader in) {
			super(in);
			this.reader = in;
		}

		@Override
		public String readLine() throws IOException {
			String line = reader.readLine();
			if(line !=null){
				buffer.append(line).append(NEW_LINE);
			}
			return line;
		}


		@Override
		public void mark(int readAheadLimit) throws IOException {
			resetPosition = buffer.length();
			reader.mark(readAheadLimit);
		}

		@Override
		public void reset() throws IOException {
			//call super first which does all the error checking
			reader.reset();
			buffer.setLength(resetPosition);
			resetPosition =-1;
		}
		
		
		public String getBufferedLines(){
			return buffer.toString();
		}
		public void resetBuffer(){
			buffer.setLength(0);
		}

		@Override
		public int read() throws IOException {
			return reader.read();
		}

		@Override
		public int read(char[] cbuf, int off, int len) throws IOException {
			return reader.read(cbuf, off, len);
		}

		@Override
		public long skip(long n) throws IOException {
			return reader.skip(n);
		}

		@Override
		public boolean ready() throws IOException {
			return reader.ready();
		}

		@Override
		public boolean markSupported() {
			return reader.markSupported();
		}

		@Override
		public void close() throws IOException {
			reader.close();
		}

		@Override
		public Stream<String> lines() {
			return reader.lines();
		}
		
	}

	@Override
	public boolean isDefault() {
		return true;
	}

}
