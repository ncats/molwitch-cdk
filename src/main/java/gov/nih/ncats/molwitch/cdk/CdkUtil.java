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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;

import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.Expr.Type;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.ChemicalSource;

public class CdkUtil {

    private static final Aromaticity AROMATICITY = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));

    public static void aromatize(IAtomContainer container) throws CDKException {
		setImplicitHydrogensIfNeeded(container, false);
//		 Aromaticity.cdkLegacy().apply(container);
        AROMATICITY.apply(container);
    }
    
    public static IAtomContainer safeClone(IAtomContainer iac) {
        //the CDK clone doesn't actually clone deep, but still allows modifications
        return (IAtomContainer)asChemical(iac).copy().getImpl().getWrappedObject();
    }
    public static Chemical asChemical(IAtomContainer iac) {
        return new Chemical(new CdkChemicalImpl(iac,(ChemicalSource)null));
    }
    
    public static void poorMansKekulize(IAtomContainer iac) {
        throw new UnsupportedOperationException("simple kekule not yet implemnted");
        
        //basic logic:
        // 1. Find all atoms with aromatic bonds
        // 2. For any atom that has only 1 aromatic bond: fail
        // 3. For any atom that has 2 aromatic bonds, look at H count
        //    and determine where double bonds could go.
        
//        for(IBond bond : iac.bonds()){
//            if(bond.isAromatic()){
//                
//            }
//        }
    }
    
    
    
	public static IAtomContainer kekulizeIfNeeded(IAtomContainer container, boolean makeCopy) throws CDKException {
		IAtomContainer possibleCopy = setImplicitHydrogensIfNeeded(container, makeCopy);
		boolean setBonds=false;
		for(IBond bond : possibleCopy.bonds()){
			if(bond.getOrder() ==null || bond.getOrder()==Order.UNSET){
				setBonds=true;
				break;
			}
		}
		if(setBonds){
			IAtomContainer ret = possibleCopy;
			if(possibleCopy ==container && makeCopy){
				try {
					ret = safeClone(container);
				} catch (Exception e) {
					//shouldn't happen container support clones
				}

			}
			Kekulization.kekulize(ret);
			for(IBond bond : ret.bonds()){
				bond.setIsAromatic(false);
			}
			return ret;
		}
		for(IBond bond : possibleCopy.bonds()){
			bond.setIsAromatic(false);
		}
		return possibleCopy;

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
	    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ret);
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
		IAtomContainer container =CdkUtil.getChemObjectBuilder().newAtomContainer();
		Smarts.parse(container, smarts);

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
		QueryAtomPerceptor.percieve(container);
		return container;
	}
	
	public static void removeQueryAtoms(IAtomContainer containerToFingerprint) {
        //iterator.remove() doesn't seem to work
        List<Integer> indexesToRemove = new ArrayList<>();
        int i=0;
        for(IAtom next : containerToFingerprint.atoms()){

//          System.out.println("atom  = " + next + " sym = " + next.getSymbol());
            if(next.getSymbol() ==null){
                indexesToRemove.add(Integer.valueOf(i));
            }
            i++;
        }
        final IAtomContainer fia=containerToFingerprint;
         //iterate in reverse to not mess up atom order of deletions
        Comparator<Integer> ci=Comparator.naturalOrder();
        indexesToRemove.stream().sorted(ci.reversed()).forEach(k->{
            fia.removeAtom(k);
        });
	}
	public static void removeQueryBonds(IAtomContainer containerToFingerprint) {

	    Iterator<IBond> bondIter = containerToFingerprint.bonds().iterator();
	    while(bondIter.hasNext()){
	        IBond next = bondIter.next();
	        System.out.println(next.getOrder());
	        if(next.getOrder() ==null || next.getOrder() == Order.UNSET){
	            if(!next.isAromatic()) {
	              if(CdkUtil.isSubtleQueryBond(next) 
	                      && !CdkUtil.simplifyQBond(next)
	                      ) {
	                    bondIter.remove();    
	              }                       
	            }
	        }
	    }
	}
	public static IAtomContainer getSimplifiedContainer(IAtomContainer iac) {
	    if(iac instanceof QueryAtomContainer) {
            QueryAtomContainer qac = (QueryAtomContainer)iac;
            CdkChemicalImpl cci = new CdkChemicalImpl(qac, (ChemicalSource)null);
            Chemical cq = new Chemical(cci);
            Chemical qq=cq.copy();
            IAtomContainer ii = (IAtomContainer)qq.getImpl().getWrappedObject();
            return ii;
        }
	    return iac;
	}


	
	private static String getSymbolForAtomExpression(Expr exp){
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
	public static void navNodes(Expr exr, BiConsumer<Integer, Expr> elist, int l){
		elist.accept(l, exr);
		if(exr.type().equals(Expr.Type.OR) || exr.type().equals(Expr.Type.AND)){
			navNodes(exr.left(), elist,l+1);
			navNodes(exr.right(), elist,l+1);
		}else if(exr.type().equals(Expr.Type.NOT)){
			navNodes(exr.left(), elist,l+1);
		}
	}
	


	public static String rep(String m, int l){
		String r="";
		for(int i=0;i<l;i++){
			r+=m;
		}
		return r;
	}

	public static boolean isPseudoAtom(IAtom atom) {
		IAtom deref = AtomRef.deref(atom);
		return deref instanceof IPseudoAtom;
	}
    public static QueryAtomContainer asQueryAtomContainer(IAtomContainer ia){
    	QueryAtomContainer qac=QueryAtomContainer.create(ia);
    	for(int i=0;i<qac.getBondCount();i++){
    		QueryBond ib=(QueryBond)qac.getBond(i);
    		
    		IBond ibo=ia.getBond(i);
    		ibo = BondRef.deref(ibo);
    		if(ibo instanceof QueryBond){
    			ib.setExpression(((QueryBond)ibo).getExpression());

    			if(ibo.isAromatic()){
    				ib.setExpression(new Expr(Expr.Type.IS_AROMATIC));
    				
    			}else if(ib.getExpression().type().equals(Expr.Type.ORDER)){
					ib.getExpression().setPrimitive(Expr.Type.ALIPHATIC_ORDER, ib.getExpression().value());
    			}
    		}else{
    			if(!ibo.isAromatic()){
    				if(ibo.getOrder()==null || ibo.getOrder().equals(Order.UNSET)){
    					ib.setExpression(new Expr(Expr.Type.TRUE));
    				}else{
    					ib.setExpression(new Expr(Expr.Type.ALIPHATIC_ORDER,ibo.getOrder().numeric()));
    				}
    			}else{
    				ib.setExpression(new Expr(Expr.Type.IS_AROMATIC));
    			}
    			
    		}
    		if(ib.getExpression().type().equals(Expr.Type.STEREOCHEMISTRY)){
    			ib.setExpression(new Expr(Expr.Type.TRUE));
    		}
    	}        	
    	for(int i=0;i<qac.getAtomCount();i++){
    		QueryAtom iat=(QueryAtom)qac.getAtom(i);
    		
    		IAtom iao=ia.getAtom(i);
    		iao = AtomRef.deref(iao);
    		if(iao instanceof QueryAtom){
    			
    			iat.setExpression(((QueryAtom)iao).getExpression());
    			iat.setSymbol(getSymbolForAtomExpression(iat.getExpression()));
    			if(iat.getExpression().type().equals(Expr.Type.ALIPHATIC_ELEMENT) ||
    					iat.getExpression().type().equals(Expr.Type.ELEMENT)||
    					iat.getExpression().type().equals(Expr.Type.AROMATIC_ELEMENT)){
    				iat.setAtomicNumber(iat.getExpression().value());
    				iat.setSymbol(PeriodicTable.getSymbol(iat.getExpression().value()));
    			}
    			
    		}else{
    			if(iao.getAtomicNumber()==null){
    				iat.setSymbol("A");
    			}else{
    				iat.setExpression(new Expr(Expr.Type.ELEMENT,iao.getAtomicNumber()));
    				iat.setSymbol(iao.getSymbol());
    			}
    		}
    		if(iao.getCharge()!=null){
    			iat.setCharge(iao.getCharge());
    		}
    		if(iao.getMassNumber()!=null){
    			iat.setMassNumber(iat.getMassNumber());
    		}

			iat.setPoint2d(iao.getPoint2d());
    	}
    	
    	return qac;
    }
    
    public static boolean isSubtleQueryBond(IBond b) {
        if(b instanceof QueryBond) {
            QueryBond qb = (QueryBond)b;
            Expr ex=qb.getExpression();
            if(ex.left()!=null || ex.right()!=null)return true;
            
            org.openscience.cdk.isomorphism.matchers.Expr.Type t=ex.type();
//            org.openscience.cdk.isomorphism.matchers.Expr.Type.
            switch(t) {
                case IS_AROMATIC:
                case ORDER:  
                    return false;
                default:
                    return true;
            }
//            ex.
//            System.out.println(t);
        }else {
            if(b.getOrder()==null || b.getOrder() == Order.UNSET) {
                return true;
            }
        }
        return false;
    }
    public static boolean simplifyQBond(IBond b) {
        if(b instanceof QueryBond) {
            QueryBond qb = (QueryBond)b;
            Expr ex=qb.getExpression();
            if(ex.left()!=null || ex.right()!=null)return false;
            switch(ex.type()) {
            // Under certain circumstances a parse of a smiles/smarts pattern
            // may result in SINGLE_OR_AROMATIC or DOUBLE_OR_AROMATIC bonds
            // when really plain single or plain double are intended
            //
            // For now, the feature to "rescue" these cases is turned off
            
//                case SINGLE_OR_AROMATIC:
//                    ex.setPrimitive(Type.ORDER, 1);
//                    return true;
//                case DOUBLE_OR_AROMATIC:
//                    ex.setPrimitive(Type.ORDER, 2);
//                    return true;
            default:
                break;
            }
        }
        return false;
    }
    
    public static boolean canWrite(IAtomContainer mol){
    	for(IAtom ia : mol.atoms()){
    		if(ia.getSymbol() == null)return false;
    	}
    	for(IBond ib : mol.bonds()){
    		if(ib.getOrder() == null || ib.getOrder().equals(Order.UNSET)) {
    		    if(!ib.isAromatic()) {
    		        return false;
    		    }
    		}
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
