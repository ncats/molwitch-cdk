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

package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.GraphInvariant;

public class CdkGraphInvariant implements GraphInvariant {

	private long[] inv;
	
	public CdkGraphInvariant(long[] inv) {
		this.inv = inv;
	}

	@Override
	public int getNumberOfAtoms() {
		return inv.length;
	}

	@Override
	public long getAtomInvariantValue(int atomIndex) {
		return inv[atomIndex];
	}

	@Override
	public boolean isAtomCompatible(int i, int j) {
		return inv[i] == inv[j];
	}

}
