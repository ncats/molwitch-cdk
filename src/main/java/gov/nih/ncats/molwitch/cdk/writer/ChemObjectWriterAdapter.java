/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2025.
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

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;
import java.util.Collection;
import java.util.function.Function;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.IChemObjectWriter;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.io.listener.IChemObjectIOListener;
import org.openscience.cdk.io.setting.IOSetting;

public class ChemObjectWriterAdapter<T extends IChemObject> implements IChemObjectWriter, Closeable{

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
