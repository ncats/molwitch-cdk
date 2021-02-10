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
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.algorithm.single.SingleMappingHandler;
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
        IAtomContainer container= CdkUtil.toAtomContainer(chemical);
        if(hasQueryAtomsOrBonds(container)){
            query = new QueryAtomContainer(container, container.getBuilder());
        }else{
            query = container;
        }
//        query = new QueryAtomContainer(container, container.getBuilder());
    }

    private static boolean hasQueryAtomsOrBonds(IAtomContainer container){
        for(IAtom a : container.atoms()){
            if(AtomRef.deref(a) instanceof IQueryAtom){
                return true;
            }
            for(IBond b : container.getConnectedBondsList(a)){
                if(b instanceof IQueryBond){
                    return true;
                }
            }
        }
        return false;
    }
    @Override
    public Optional<int[]> search(Chemical targetChemical) {
        //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
  //Algorithm is VF2MCS
  //Bond Sensitive is set True
  //Ring Match is set True
        IAtomContainer target = CdkUtil.toAtomContainer(targetChemical);
        try {
            //set matchBonds to false if set to true a substructure with extra bonds fails test ?
            Substructure smsd = new Substructure(query, target, false, false, true, false, false);

            if(!smsd.isSubgraph()){
                return Optional.empty();
            }
            AtomAtomMapping atomMapping;
            if(query.getAtomCount() ==1 || target.getAtomCount()==1){
                //for some reason SDSM doesn't report mapping
                //of single atom but it does set isSubgraph if they match

                SingleMappingHandler mcs;
                if (!(query instanceof IQueryAtomContainer) && !(target instanceof IQueryAtomContainer)) {
                    mcs = new SingleMappingHandler(query, target, false);
                } else {
                    mcs = new SingleMappingHandler((IQueryAtomContainer) query,target);
                }
                atomMapping = mcs.getFirstAtomMapping();
//                    return mcs.getAllAtomMapping() != null && !mcs.getAllAtomMapping().isEmpty();

            }else {
                atomMapping = smsd.getFirstAtomMapping();
            }
            if(atomMapping.isEmpty()){
                return Optional.empty();
            }
            int[] map = new int[query.getAtomCount()];
            for (Map.Entry<IAtom, IAtom> mapping : atomMapping.getMappingsByAtoms().entrySet()) {
                IAtom sourceAtom = mapping.getKey();
                IAtom targetAtom = mapping.getValue();
                int queryMappingNumber = atomMapping.getQueryIndex(sourceAtom);
                //Get the mapped atom number in Target AtomContainer
                int targetMappingNumber = atomMapping.getTargetIndex(targetAtom);
                map[queryMappingNumber] = targetMappingNumber;
//           System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
//           System.out.println(queryMappingNumber + " " + targetMappingNumber);
            }
            return Optional.ofNullable(map);
        } catch (CDKException e) {
            e.printStackTrace();
            return Optional.empty();
        }
//        Isomorphism comparison = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);


    }
}
