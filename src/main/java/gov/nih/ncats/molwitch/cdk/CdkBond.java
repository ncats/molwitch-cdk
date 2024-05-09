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

package gov.nih.ncats.molwitch.cdk;


import java.util.Optional;

import org.openscience.cdk.BondRef;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

import gov.nih.ncats.common.sneak.Sneak;
import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;

public class CdkBond implements Bond{

	private final IBond bond;
	private final CdkChemicalImpl parent;


	public CdkBond(IBond bond, CdkChemicalImpl parent) {
		this.bond = bond;
		this.parent = parent;
	}
	
	




	@Override
	public void setBondType(BondType type) {
		switch(type){
		case SINGLE : bond.setOrder(Order.SINGLE);
						break;
		case DOUBLE  : bond.setOrder(Order.DOUBLE);
		break;
		case TRIPLE : bond.setOrder(Order.TRIPLE);
		break;
		case QUADRUPLE : bond.setOrder(Order.QUADRUPLE);
						break;
		//CDK has aromatic as a flag not a bond...
						
		case AROMATIC : bond.setIsAromatic(true);
						break;
		}
		
	}





	@Override
	public boolean isQueryBond() {
		return BondRef.deref(bond) instanceof IQueryBond;
	}



	@Override
	public Bond switchParity() {
		IBond clone = null;
		try {
			clone = (IBond) bond.clone();
		} catch (CloneNotSupportedException e) {
			Sneak.sneakyThrow(e);
		}
		
		clone.setOrder(toOrder(getBondType().switchParity()));
		return new CdkBond(clone, parent);
	}



	@Override
	public Atom getOtherAtom(Atom a) {
		IAtom atom =((CdkAtom)a).getAtom();
		IAtom other = bond.getOther(atom);
		return parent.getCdkAtomFor(other);
	}

	@Override
	public Atom getAtom1() {
		//is this zero based? or 1 based?
		IAtom atom = bond.getAtom(0);		
		return parent.getCdkAtomFor(atom);
	}

	@Override
	public Atom getAtom2() {
		//is this zero based? or 1 based?
		IAtom atom = bond.getAtom(1);
		return parent.getCdkAtomFor(atom);
	}

	@Override
	public BondType getBondType() {
		
		if(bond.isAromatic()){
			if(parent.isAromatic()){
				return BondType.AROMATIC;
			}
		}
		if(bond.getOrder() ==null) {
		    return null;
		}
		switch(bond.getOrder()){
			case SINGLE : return BondType.SINGLE;
			case DOUBLE : return BondType.DOUBLE;
			case TRIPLE : return BondType.TRIPLE;
			case QUADRUPLE : return BondType.QUADRUPLE;
			case UNSET : {
				//sometimes CDK marks an aromatic bond as unset
				if(bond.getAtom(0).isAromatic() 
						&& bond.getAtom(1).isAromatic() ) {
					return BondType.AROMATIC;
				}
			}
			default:{
					
					return null;
			}
		}
	
	}
	
	private Order toOrder(BondType type){
		switch(type){
			case SINGLE : return Order.SINGLE;
			case DOUBLE : return Order.DOUBLE;
			case TRIPLE : return Order.TRIPLE;
			case QUADRUPLE : return Order.QUADRUPLE;
			default : return Order.DOUBLE; //TODO what do we do here?
		}
	}

	IBond getBond() {
		return bond;
	}
	

	@Override
	public Stereo getStereo() {
		IBond.Stereo cdkStereo = bond.getStereo();
		if(cdkStereo == null){
			return null;
		}
		switch(cdkStereo){
			case DOWN : return Stereo.DOWN;
			case DOWN_INVERTED : return Stereo.DOWN_INVERTED;
			case NONE : return Stereo.NONE;
			case UP : return Stereo.UP;
			case UP_INVERTED : return Stereo.UP_INVERTED;
			case UP_OR_DOWN : return Stereo.UP_OR_DOWN;
			case UP_OR_DOWN_INVERTED : return Stereo.UP_OR_DOWN_INVERTED;

			default :
				return Stereo.NONE;
		}
	}

	@Override
	public void setStereo(Stereo stereo) {
		if(stereo ==null){
			bond.setStereo(null);
			return;
		}
		switch (stereo) {
		case DOWN:
			bond.setStereo(IBond.Stereo.DOWN);
			break;
		case DOWN_INVERTED:
			bond.setStereo(IBond.Stereo.DOWN_INVERTED);
			break;
		case NONE:
			bond.setStereo(IBond.Stereo.NONE);
			break;
		case UP:
			bond.setStereo(IBond.Stereo.UP);
			break;
		case UP_INVERTED:
			bond.setStereo(IBond.Stereo.UP_INVERTED);
			break;
		case UP_OR_DOWN:
			bond.setStereo(IBond.Stereo.UP_OR_DOWN);
			break;
		case UP_OR_DOWN_INVERTED:
			bond.setStereo(IBond.Stereo.UP_OR_DOWN_INVERTED);
			break;
		}
		parent.setDirty();
		
	}

	@Override
	public DoubleBondStereo getDoubleBondStereo() {
		if(getBondType() != BondType.DOUBLE){
			return null;
		}
		parent.cahnIngoldPrelogSupplier.get();
		String value = Optional.ofNullable(bond.getProperty(CDKConstants.CIP_DESCRIPTOR)).map(bt->bt.toString()).orElse(null);
		if("Z".equals(value)) {
			 return DoubleBondStereo.Z_CIS;
		}
		if("E".equals(value)) {
			 return DoubleBondStereo.E_TRANS;
		}
		if(bond.getStereo()== IBond.Stereo.E_OR_Z){
			return DoubleBondStereo.E_OR_Z;
		}
		return DoubleBondStereo.NONE;
		
		
		
		
	}



	@Override
	public boolean isAromatic() {
		return bond.isAromatic();
	}






	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;

		CdkBond cdkBond = (CdkBond) o;

		if (!bond.equals(cdkBond.bond)) return false;
		return parent.equals(cdkBond.parent);

	}

	@Override
	public int hashCode() {
		int result = bond.hashCode();
		result = 31 * result + parent.hashCode();
		return result;
	}



	@Override
	public boolean isInRing() {
		try {
			parent.ringsSearcherSupplier.get();
			return bond.isInRing();
		}
		catch(Exception ex){
			System.err.println("Error in isInRing()");
			ex.printStackTrace();
		}
		return false;
	}


    @Override
    public String toString() {
        return "CdkBond{" +
                "bond=" + bond +
                '}';
    }

    public static IBond getIBondFor(Bond b) {
		return ((CdkBond)b).bond;
		
	}
}
