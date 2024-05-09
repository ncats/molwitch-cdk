/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2024.
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

package gov.nih.ncats.molwitch.cdk.fingerprinters;

import java.util.Collections;
import java.util.Set;

import org.openscience.cdk.fingerprint.ExtendedFingerprinter;

import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.FingerprintSpecification;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.PathBasedSpecification;
import gov.nih.ncats.molwitch.spi.FingerprinterImpl;

public class DefaultFingerPrinter implements FingerprinterImpl{

	@Override
	public boolean supports(FingerprintSpecification options) {
		return Fingerprinters.FingerprintSpecification.PATH_BASED.name().equals(options.name());
	}

	@Override
	public boolean isDefault() {
		return true;
	}

	@Override
	public Set<String> getSupportedAlgorithmNames() {
		return Collections.singleton(Fingerprinters.FingerprintSpecification.PATH_BASED.name());
	}

	@Override
	public Fingerprinter createDefaultFingerprinter() {
		return createFingerPrinterFor(null);
	}

	@Override
	public Fingerprinter createFingerPrinterFor(FingerprintSpecification fingerPrinterOptions) {
		FingerprinterAdapter adapter;
		if(fingerPrinterOptions instanceof PathBasedSpecification) {
		PathBasedSpecification options = (PathBasedSpecification)fingerPrinterOptions;

			adapter = new FingerprinterAdapter( new org.openscience.cdk.fingerprint.Fingerprinter(options.getLength(), options.getDepth()));
		}else{
			adapter = new FingerprinterAdapter( new org.openscience.cdk.fingerprint.Fingerprinter());
		}
		adapter.setRemoveQueryAtomsAndBonds(true);
		//We don't want this true by default
//		adapter.setForceExplicitH(true);
		return adapter;
	}

}
