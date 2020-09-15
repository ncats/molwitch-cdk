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
import gov.nih.ncats.molwitch.cdk.CdkAtom;
import gov.nih.ncats.molwitch.cdk.CdkUtil;
import gov.nih.ncats.molwitch.search.MolSearcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;

import java.io.IOException;
import java.util.Map;
import java.util.Optional;

public class CdkMolSearcher implements MolSearcher {
    IAtomContainer query;
    public CdkMolSearcher(String smartsQuery){
        try {
            query = CdkUtil.parseSmarts(smartsQuery);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
    public CdkMolSearcher(Chemical chemical){
        query = CdkUtil.toAtomContainer(chemical);
    }
    @Override
    public Optional<int[]> search(Chemical targetChemical) {
        //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
  //Algorithm is VF2MCS
  //Bond Sensitive is set True
  //Ring Match is set True
        IAtomContainer target = CdkUtil.toAtomContainer(targetChemical);
  Isomorphism comparison = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);
            int[] map = new int[query.getAtomCount()];


           for (Map.Entry<IAtom, IAtom> mapping : comparison.getFirstAtomMapping().getMappingsByAtoms().entrySet()) {
               IAtom sourceAtom = mapping.getKey();
               IAtom targetAtom = mapping.getValue();
               map[sourceAtom.getIndex()] = targetAtom.getIndex();
    //           System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
    //           System.out.println(atomatomMapping.getQueryIndex(sourceAtom) + " " + atomatomMapping.getTargetIndex(targetAtom));
           }
        return Optional.ofNullable(map);
    }
}
