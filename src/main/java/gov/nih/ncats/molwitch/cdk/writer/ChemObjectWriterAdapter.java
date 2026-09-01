package gov.nih.ncats.molwitch.cdk.writer;

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;
import java.util.Collection;
import java.util.function.Function;
import java.util.logging.Logger;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.IChemObjectWriter;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.io.listener.IChemObjectIOListener;
import org.openscience.cdk.io.setting.IOSetting;

public class ChemObjectWriterAdapter<T extends IChemObject> implements IChemObjectWriter, Closeable{
	private static Logger logger = Logger.getLogger("ChemObjectWriterAdapter");

	private final IChemObjectWriter delegate;
	private final Function<T, T> adapter;
	public static <T extends IChemObject> ChemObjectWriterAdapter<T> create(IChemObjectWriter delegate, Function<T, T> adapter){
		return new ChemObjectWriterAdapter<>(delegate, adapter);
	}
	private ChemObjectWriterAdapter(IChemObjectWriter delegate, Function<T, T> adapter) {
		this.delegate = delegate;
		this.adapter = adapter;
	}

	@Override
	public IResourceFormat getFormat() {
		return delegate.getFormat();
	}

	@Override
	public boolean accepts(Class<? extends IChemObject> classObject) {
		return delegate.accepts(classObject);
	}

	@Override
	public void close() throws IOException {
		delegate.close();
	}

	@Override
	public IOSetting[] getIOSettings() {
		return delegate.getIOSettings();
	}

	@Override
	public void addChemObjectIOListener(IChemObjectIOListener listener) {
		delegate.addChemObjectIOListener(listener);
	}

	@Override
	public void removeChemObjectIOListener(IChemObjectIOListener listener) {
		delegate.removeChemObjectIOListener(listener);
	}

	@Override
	public Collection<IChemObjectIOListener> getListeners() {
		return delegate.getListeners();
	}

	@Override
	public <S extends IOSetting> S addSetting(IOSetting setting) {
		return delegate.addSetting(setting);
	}

	@Override
	public void addSettings(Collection<IOSetting> settings) {
		delegate.addSettings(settings);
	}

	@Override
	public boolean hasSetting(String name) {
		return delegate.hasSetting(name);
	}

	@Override
	public <S extends IOSetting> S getSetting(String name) {
		return delegate.getSetting(name);
	}

	@Override
	public <S extends IOSetting> S getSetting(String name, Class<S> c) {
		return delegate.getSetting(name, c);
	}

	@Override
	public Collection<IOSetting> getSettings() {
		return delegate.getSettings();
	}

	@SuppressWarnings("unchecked")
	@Override
	public void write(IChemObject object) throws CDKException {
		if( delegate instanceof MDLV2000Writer) {
			logger.info("in write, will cast");
			((MDLV2000Writer) delegate).customizeJob();
		} else if( delegate instanceof SDFWriter) {
			logger.info("in write, will cast as SDFWriter");
			((SDFWriter) delegate).customizeJob();
		} else {
			logger.info("in write, cast NOT valid. Object is a " + delegate.getClass().getName());
		}
		delegate.write(adapter.apply((T)object));
		
	}

	@Override
	public void setWriter(Writer writer) throws CDKException {
		delegate.setWriter(writer);
	}

	@Override
	public void setWriter(OutputStream writer) throws CDKException {
		delegate.setWriter(writer);
	}

}
