package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.ChemicalBuilder;
import gov.nih.ncats.molwitch.ChemicalSource.Type;
import gov.nih.ncats.molwitch.inchi.InChiResult;
import gov.nih.ncats.molwitch.inchi.InChiResult.Status;
import gov.nih.ncats.molwitch.internal.source.StringSource;
import gov.nih.ncats.molwitch.spi.InchiImplFactory;
import io.github.dan2097.jnainchi.InchiOptions;
import io.github.dan2097.jnainchi.InchiStatus;
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;

import java.io.IOException;

public class CdkChemicalInchiImplFactory implements InchiImplFactory{

	private final InChIGeneratorFactory factory;
	//CDK by default turns of auxInfo
	private static InchiOptions MOLWITCH_INCHI_OPTIONS = new InchiOptions.InchiOptionsBuilder().build();
	public CdkChemicalInchiImplFactory() {
		try {
			factory = InChIGeneratorFactory.getInstance();
		} catch (CDKException e) {
			throw new IllegalStateException("could not initialize Inchi generator", e);
		}
		
	}
	private Chemical handleQueryAtoms(Chemical m) throws IOException {

		if(!m.atoms().filter(CdkChemicalInchiImplFactory::isAtomThatMessesUpInchi).findAny().isPresent()){
			//no query atoms do nothing
			return m;
		}
		Chemical copy = m.copy();
		for (Atom a : copy.getAtoms()) {
			if (isAtomThatMessesUpInchi(a)) {
				// this is what marvinjs specifies for atom *
				a.setAtomicNumber(2); // force this to be helium

			}
		}

		return copy;
	}

	private static boolean isAtomThatMessesUpInchi(Atom atom) {
		IAtom a = AtomRef.deref(CdkAtom.getIAtomFor(atom));
		boolean ret =  a instanceof IPseudoAtom || a instanceof IQueryAtom;
		return ret;
	}

	@Override
	public InChiResult asStdInchi(Chemical chemical, boolean trustCoordinates) throws IOException {
		try {

			Chemical ichem=handleQueryAtoms(chemical);
			ichem.kekulize();
			
			InChIGenerator gen = factory.getInChIGenerator(CdkUtil.toAtomContainer(ichem), MOLWITCH_INCHI_OPTIONS);
		
			InChiResult.Status status = toChemkitStatus(gen.getStatus());
			String inchi = gen.getInchi();
			if(inchi ==null){
				return new InChiResult.Builder(status)
						.setAuxInfo(gen.getAuxInfo()==null?"":gen.getAuxInfo())
						.setInchi("")
						.setKey("")
						.setMessage(gen.getLog()==null? "":gen.getLog())
						.build();
			}
			return new InChiResult.Builder(status)
						.setAuxInfo(gen.getAuxInfo()==null?"":gen.getAuxInfo())
						.setInchi(gen.getInchi()==null?"":gen.getInchi())
						.setKey(gen.getInchiKey()==null?"" :gen.getInchiKey())
						.setMessage(gen.getLog()==null? "":gen.getLog())
						.build();
			
		
		} catch (Throwable e) {
			e.printStackTrace();
			throw new IOException("error computing Inchi for " + chemical.toSmarts(), e);
		} 
	}

	private InChiResult.Status toChemkitStatus(InchiStatus returnStatus) {
		if(returnStatus == InchiStatus.SUCCESS) {
			return Status.VALID;
		}
		if(returnStatus == InchiStatus.WARNING) {
			return Status.WARNING;
		}
		return Status.ERROR;
	}
	@Override
	public Chemical parseInchi(String inchi) throws IOException {
		try {
			InChIToStructure toStruc= factory.getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());
		
			if(toStruc.getStatus() == InchiStatus.SUCCESS || toStruc.getStatus() == InchiStatus.WARNING) {
				
				return ChemicalBuilder._fromImpl(new CdkChemicalImpl( toStruc.getAtomContainer(), new StringSource(inchi, Type.INCHI)))
								.computeCoordinates(true)
//								.computeStereo(true)
								.build();
							
			}
			throw new IOException("error trying to parse inchi '"+ inchi + "' : " + toStruc.getMessage());
		
		} catch (CDKException e) {
			throw new IOException("error trying to parse inchi '"+ inchi + "'", e);
		}
	}

}
