/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2023.
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

import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryBond;

public class QueryAtomPerceptor {

	
	public static void percieve(IAtomContainer container) throws IOException{
		for(IAtom atom : container.atoms()){
			IAtom deref = AtomRef.deref(atom);
			if(deref instanceof QueryAtom ){
				Expr expr= ((QueryAtom)deref).getExpression();
//				System.out.println(expr.type());
				if(expr.type() == Expr.Type.ELEMENT || expr.type() == Expr.Type.ALIPHATIC_ELEMENT || expr.type() == Expr.Type.AROMATIC_ELEMENT){
					int atomicNumber = expr.value();
					IElement element = Isotopes.getInstance().getElement(atomicNumber);
					atom.setAtomicNumber(atomicNumber);
					atom.setSymbol(element.getSymbol());
					if(expr.type() == Expr.Type.AROMATIC_ELEMENT){
						atom.setIsAromatic(true);
					}
				}
			}
		}
		
		for(IBond bond : container.bonds()){
			IBond deref = BondRef.deref(bond);
//			System.out.println(deref);
			if(deref instanceof QueryBond){
				Expr expr = ((QueryBond)deref).getExpression();
//				System.out.println("\t" + expr);
				if(expr.type() == Expr.Type.ALIPHATIC_ORDER || expr.type() == Expr.Type.ORDER){
					switch(expr.value()){
						case 1: bond.setOrder(Order.SINGLE);
								break;
						case 2: bond.setOrder(Order.DOUBLE);
								break;
						case 3: bond.setOrder(Order.TRIPLE);
							break;
						case 4: bond.setOrder(Order.QUADRUPLE);
							break;
					}
					
				}
				if(expr.type() == Expr.Type.SINGLE_OR_AROMATIC){
					//this is probably an aromatic ring?
					bond.setOrder(Order.SINGLE);
					bond.setIsAromatic(true);
					bond.setIsInRing(true);
				}
			}
		}
	}
}
