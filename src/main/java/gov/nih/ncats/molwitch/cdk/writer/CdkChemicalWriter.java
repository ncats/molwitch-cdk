package gov.nih.ncats.molwitch.cdk.writer;

import java.io.IOException;
import java.util.Properties;
import java.util.logging.Logger;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.IChemObjectWriter;

import gov.nih.ncats.molwitch.cdk.CdkChemicalImpl;
import gov.nih.ncats.molwitch.cdk.CdkUtil;
import gov.nih.ncats.molwitch.spi.ChemicalImpl;
import gov.nih.ncats.molwitch.spi.ChemicalWriterImpl;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.listener.PropertiesListener;

public class CdkChemicalWriter implements ChemicalWriterImpl{
	private static Logger logger = Logger.getLogger("CdkChemicalWriter");

	private final IChemObjectWriter writer;

	public CdkChemicalWriter(IChemObjectWriter writer) {
		this.writer = writer;
	}

	@Override
	public void close() throws IOException {
		writer.close();
	}
    
	@Override
	public void write(ChemicalImpl impl) throws IOException {
		CdkChemicalImpl chem =(CdkChemicalImpl)impl;
		IAtomContainer mol =CdkUtil.getUsableFormOfAtomContainer(chem.getContainer());
		Properties sdfWriterProps = new Properties();
		sdfWriterProps.put("WriteAromaticBondTypes", "true");
		writer.addChemObjectIOListener(new PropertiesListener(sdfWriterProps));
		try {
			writer.write(mol);
		}catch(Throwable e) {
			throw new IOException("error writing container " + mol.getID(), e);
			
		}
		
	}

	
	
}
