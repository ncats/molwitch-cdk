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

package gov.nih.ncats.molwitch.cdk.fingerprinters;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.openscience.cdk.fingerprint.CircularFingerprinter;

import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.AbstractCfpOptions;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.EcfpSpecification;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.FingerprintSpecification;
import gov.nih.ncats.molwitch.spi.FingerprinterImpl;

public class EcfpFingerPrinter implements FingerprinterImpl{

	private static List<String> supportedAlgorithms = Arrays.asList(Fingerprinters.FingerprintSpecification.ECFP.name(),
										Fingerprinters.FingerprintSpecification.FCFP.name()); 
	public EcfpFingerPrinter() {
	}

	@Override
	public boolean supports(FingerprintSpecification options) {
		return supportedAlgorithms.contains(options.name());
	}

	@Override
	public boolean isDefault() {
		return false;
	}

	@Override
	public Set<String> getSupportedAlgorithmNames() {
		return Collections.unmodifiableSet(new HashSet<>(supportedAlgorithms));
	}

	@Override
	public Fingerprinter createFingerPrinterFor(FingerprintSpecification fingerPrinterOptions) {
		AbstractCfpOptions<?> ecfpOptions = (AbstractCfpOptions<?>)fingerPrinterOptions;
		int diameter = ecfpOptions.getDiameter();
		int classType;
		if(ecfpOptions instanceof EcfpSpecification) {
			switch(diameter){
				case 0 : classType= CircularFingerprinter.CLASS_ECFP0;
						break;
				case 2:  classType= CircularFingerprinter.CLASS_ECFP2;
						break;
				case 4:  classType= CircularFingerprinter.CLASS_ECFP4;
				break;
				default: classType= CircularFingerprinter.CLASS_ECFP6;
			}
		}else {
			switch(diameter){
			case 0 : classType= CircularFingerprinter.CLASS_FCFP0;
					break;
			case 2:  classType= CircularFingerprinter.CLASS_FCFP2;
					break;
			case 4:  classType= CircularFingerprinter.CLASS_FCFP4;
			break;
			default: classType= CircularFingerprinter.CLASS_FCFP6;
		}
		}
		return new FingerprinterAdapter(new CircularFingerprinter(classType, ecfpOptions.getBitLength()));
	}
}
