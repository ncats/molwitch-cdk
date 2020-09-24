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

package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Chemical;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CdkUtil {

    private static final Aromaticity AROMATICITY = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));

    public static void aromatize(IAtomContainer container) throws CDKException {
		setImplicitHydrogensIfNeeded(container, false);
        AROMATICITY.apply(container);

    }

    public static IAtomContainer setImplicitHydrogensIfNeeded(IAtomContainer container, boolean makeCopy) throws CDKException {
    	boolean recompute=false;
    	for(IAtom atom : container.atoms()){
    		if(atom.getImplicitHydrogenCount() ==null){
    			recompute=true;
    			break;
			}
		}
    	if(!recompute){
    		return container;
		}
    	IAtomContainer ret = container;
    	if(makeCopy){
			try {
				ret = container.clone();
			} catch (CloneNotSupportedException e) {
				//shouldn't happen container support clones
			}
		}
		CDKHydrogenAdder.getInstance(ret.getBuilder()).addImplicitHydrogens(ret);
    	return ret;
	}
    public static IAtomContainer toAtomContainer(Chemical chemical){
		return (IAtomContainer) chemical.getImpl().getWrappedObject();
	}
	public static IChemObjectBuilder getChemObjectBuilder() {
		return SilentChemObjectBuilder.getInstance();
	}

	public static IAtomContainer parseSmarts(String smarts) throws CDKException, IOException {
		IAtomContainer container = CdkUtil.getChemObjectBuilder().newAtomContainer();
		Smarts.parse(container, smarts);

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
		QueryAtomPerceptor.percieve(container);
		return container;
	}
	

	
	public static String getSymbolForAtomExpression(Expr exp){
		List<Expr> elist = new ArrayList<>();
		getLeafNodes(exp, elist);
		if(elist.size()>0 && elist.stream()
			 .allMatch(ex->ex.type().equals(Expr.Type.ELEMENT))){
			return "L";
		}else{
			return "A";
		}
	}
	
	public static void getLeafNodes(Expr exr, List<Expr> elist){
		if(exr.type().equals(Expr.Type.OR) || exr.type().equals(Expr.Type.AND)){
			getLeafNodes(exr.left(), elist);
			getLeafNodes(exr.right(), elist);
		}else if(exr.type().equals(Expr.Type.NOT)){
			getLeafNodes(exr.left(), elist);
		}else{
			elist.add(exr);
		}
	}
	

    
    public static QueryAtomContainer asQueryAtomContainer(IAtomContainer ia){
    	QueryAtomContainer qac=QueryAtomContainer.create(ia);
    	for(int i=0;i<qac.getBondCount();i++){
    		QueryBond ib=(QueryBond)qac.getBond(i);
    		
    		IBond ibo=ia.getBond(i);
    		if(ibo instanceof QueryBond){
    			ib.setExpression(((QueryBond)ibo).getExpression());
    		}else{
    			if(!ibo.isAromatic()){
    				ib.setExpression(new Expr(Expr.Type.ALIPHATIC_ORDER,ibo.getOrder().numeric()));
    			}else{
    				ib.setExpression(new Expr(Expr.Type.IS_AROMATIC));
    			}
    		}
    	}        	
    	for(int i=0;i<qac.getAtomCount();i++){
    		QueryAtom iat=(QueryAtom)qac.getAtom(i);
    		
    		IAtom iao=ia.getAtom(i);
    		if(iao instanceof QueryAtom){
    			
    			iat.setExpression(((QueryAtom)iao).getExpression());
    			iat.setSymbol(getSymbolForAtomExpression(iat.getExpression()));
    		}else{
    			iat.setExpression(new Expr(Expr.Type.ELEMENT,iao.getAtomicNumber()));
    			iat.setSymbol(iao.getSymbol());
    		}

			iat.setPoint2d(iao.getPoint2d());
    	}
    	
    	
    	return qac;
    }
    
    public static boolean canWrite(IAtomContainer mol){
    	for(IAtom ia : mol.atoms()){
    		if(ia.getSymbol() == null)return false;
    	}
    	return true;
    }
    
    public static IAtomContainer getUsableFormOfAtomContainer(IAtomContainer mol){
    	if(!canWrite(mol)){
    		return asQueryAtomContainer(mol);
    	}
    	return mol;
    }
}
