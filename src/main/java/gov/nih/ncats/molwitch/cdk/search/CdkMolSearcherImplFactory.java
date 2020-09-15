/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2020.
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

package gov.nih.ncats.molwitch.cdk.search;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.search.MolSearcher;
import gov.nih.ncats.molwitch.spi.MolSearcherImplFactory;

public class CdkMolSearcherImplFactory implements MolSearcherImplFactory {
    @Override
    public MolSearcher create(Chemical query) {
        return new CdkMolSearcher(query);
    }

    @Override
    public MolSearcher create(String smartsPattern) {
        return new CdkMolSearcher(smartsPattern);
    }
}
