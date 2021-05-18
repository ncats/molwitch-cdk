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
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.function.Supplier;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.openscience.cdk.AtomRef;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.SingleElectron;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.RGroupQuery;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.AtomCoordinates;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chirality;
import uk.ac.ebi.beam.Element;

public class CdkAtom implements Atom{

	private static IsotopeFactory isotopeFactory;
	
	static {
		try {
		isotopeFactory = Isotopes.getInstance();
		}catch(IOException e) {
			throw new IllegalStateException("error loading isotope data", e);
		}
	}
	
	private IAtom atom;
	private CdkChemicalImpl parent;
	
	public static IAtom getIAtomFor(Atom a){
		return ((CdkAtom)a).atom;
	}
	
	public CdkAtom(IAtom atom, CdkChemicalImpl parent) {
		Objects.requireNonNull(atom);
		
		this.parent = parent;
		this.atom = atom;
		
	}
	
	

	@Override
	public boolean isValidAtomicSymbol() {
		return Element.ofSymbol(getSymbol()) !=null;
	}

	@Override
	public int getSmallestRingSize() {
		//this should set the ISINRING flags
		parent.ringsSearcherSupplier.get();
		if( !atom.isInRing()) {
			return 0;
		};
	
		
		IRingSet ringSet = Cycles.sssr(((IAtomContainer) parent.getWrappedObject())).toRingSet();
		if(!ringSet.contains(atom)) {
			return 0;
		}
		int min = Integer.MAX_VALUE;
		for(IAtomContainer r : ringSet.atomContainers()) {
			if(r.contains(atom)) {
				int len = r.getAtomCount();
				if(len < min) {
					min = len;
				}
			}
		}
		
		return min;
	}



	@Override
	public boolean isQueryAtom() {
		return AtomRef.deref(atom) instanceof IQueryAtom;
	}



	@Override
	public OptionalInt getAtomToAtomMap() {
		Object atomAtomMapping = atom.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
		if(atomAtomMapping ==null) {
			return OptionalInt.empty();
		}
		int value=0;
		if (atomAtomMapping instanceof String) {
            	//this will throw a IllegalArgumentException (runtime)
            value = Integer.parseInt((String) atomAtomMapping);
		}else if(atomAtomMapping instanceof Integer) {
			value = (Integer) atomAtomMapping;
		}
		if(value ==0) {
			return OptionalInt.empty();
		}
		return OptionalInt.of(value);
	}

	@Override
	public void setAtomToAtomMap(int value) {
		if(value ==0) {
			atom.removeProperty(CDKConstants.ATOM_ATOM_MAPPING);
		}else {
			atom.setProperty(CDKConstants.ATOM_ATOM_MAPPING, value);
		}
		
	}

	@Override
	public boolean isInRing() {
		//this should set the ISINRING flags
		parent.ringsSearcherSupplier.get();
		return atom.isInRing();
	}

	@Override
	public String getSymbol() {
		return atom.getSymbol();
	}
	
	

	@Override
	public boolean isIsotope() {
		Integer mass = atom.getMassNumber();
		if(mass ==null) {
			return false;
		}

		double majorMass = isotopeFactory.getMajorIsotopeMass(atom.getAtomicNumber());
		if(majorMass == 2D* getAtomicNumber()){
			//no major mass
			return true;
		}
		return !(Double.valueOf(mass).equals(majorMass));
	}

	@Override
	public List<? extends Bond> getBonds() {
		return parent.getBondsFor(atom);
	}

	IAtom getAtom() {
		return atom;
	}

	@Override
	public int getAtomicNumber() {
		return atom.getAtomicNumber();
	}

	@Override
	public void setAtomicNumber(int atomicNumber) {
		atom.setAtomicNumber(atomicNumber);
	}

	@Override
	public boolean hasAromaticBond() {
		return atom.isAromatic();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + atom.hashCode();
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		CdkAtom other = (CdkAtom) obj;
		if (!atom.equals(other.atom))
			return false;
		return true;
	}

	@Override
	public int getCharge() {
		Integer charge = atom.getFormalCharge();
		//assume unset = 0 ?
		if(charge ==null || CDKConstants.UNSET == charge){
			return 0;
		}
		return charge.intValue();
	}

	@Override
	public int getRadical() {
		List<ISingleElectron> list = parent.getContainer().getConnectedSingleElectronsList(atom);
		if(list.isEmpty()) {
			return 0;
		}
		return list.stream().map(ise -> ise.getElectronCount()).filter(Objects::nonNull).mapToInt(Integer::intValue).sum();
	}

	@Override
	public void setRadical(int radical) {

		//check old radical value don't do anything if match
		int oldValue = getRadical();
		if(oldValue == radical){
			return;
		}
		ISingleElectron ise = new SingleElectron(atom);
		ise.setElectronCount(radical);

		IAtomContainer container = parent.getContainer();
		container.getConnectedSingleElectronsList(atom)
								.forEach(container::removeSingleElectron);

		container.addSingleElectron(ise);

	}


	@Override
	public AtomCoordinates getAtomCoordinates() {
		Point3d points = atom.getPoint3d();
		if(points !=null) {
			return AtomCoordinates.valueOf(points.x, points.y, points.z);
		}
		Point2d otherPoints = atom.getPoint2d();
		if(otherPoints !=null) {
			return AtomCoordinates.valueOf(otherPoints.x, otherPoints.y);
		}
		return null;
	}

	@Override
	public void setAtomCoordinates(AtomCoordinates atomCoordinates) {
		if(atomCoordinates ==null) {
			atom.setPoint2d(null);
			atom.setPoint3d(null);
		}else if(atomCoordinates.is3D()) {
			atom.setPoint2d(null);
			atom.setPoint3d(new Point3d(atomCoordinates.getX(), atomCoordinates.getY(), atomCoordinates.getZ().getAsDouble()));
		}else {
			atom.setPoint3d(null);
			atom.setPoint2d(new Point2d(atomCoordinates.getX(), atomCoordinates.getY()));
		
		}
		
	}

	
	

	@Override
	public Chirality getChirality() {
		parent.cahnIngoldPrelogSupplier.get();
		String value = atom.getProperty(CDKConstants.CIP_DESCRIPTOR);
		if("R".equals(value)) {
			return Chirality.R;
		}
		if("S".equals(value)) {
			return Chirality.S;
		}
		if("EITHER".equals(value)) {
			return Chirality.Parity_Either;
		}
		
		return Chirality.Non_Chiral;
	}

	@Override
	public void setChirality(Chirality chirality) {

	}

	@Override
	public int getMassNumber() {
		Integer num= atom.getMassNumber();
		if(num==null){
			return 0;
		}
		return num.intValue();
	}
	

	@Override
	public double getExactMass() {
		// TODO Auto-generated method stub
		Double d= atom.getExactMass();
		if(d ==null) {
			return 0;
		}
		return d.doubleValue();
	}

	@Override
	public void setCharge(int charge) {
		atom.setFormalCharge(charge);

//		recomputeImplicitHydrogens();
	}

	private void recomputeImplicitHydrogens() {
		atom.setImplicitHydrogenCount(null);
		parent.setImplicitHydrogens(atom);
	}

	@Override
	public void setMassNumber(int mass) {
		if(mass <0){
			throw new IllegalArgumentException("mass can not be negative");
		}
		if(mass ==0){
			atom.setMassNumber(null);
		}else{
			atom.setMassNumber(mass);
		}
//		recomputeImplicitHydrogens();
	}

	@Override
	public int getImplicitHCount() {
		Integer ret = atom.getImplicitHydrogenCount();
		if(ret ==null){
			//bug in hydrogen count
			if(isQueryAtom()){
				return 0;
			}
			for(Bond b : getBonds()){
				if(b.isQueryBond()){
					return 0;
				}
			}
			//compute it
			return parent.setImplicitHydrogens(atom);
			
		}
		return ret;
	}
	

	@Override
	public void setImplicitHCount(Integer implicitH) {
		atom.setImplicitHydrogenCount(implicitH);
	}



	@Override
	public OptionalInt getValence() {
		parent.perceiveAtomTypesOfNonQueryAtoms.get();
		
		Integer valence= atom.getValency();
		if(valence ==null) {
			return OptionalInt.empty();
		}
		return OptionalInt.of(valence);
	}
	

	@Override
	public boolean hasValenceError() {
		//after perceiving atom types valence should not be null unless there's an error?
		//TODO maybe not set on query atoms?
		OptionalInt opin=getValence();
		if(!opin.isPresent())return true;
		
		int v = opin.getAsInt();
		
		// TODO: Need a real valence model check here ...
		// CDK probably has one. Generally specific elements
		// are only allowed to have certain valences. That is,
		// there are only certain bond count/charge/radical
		// configurations that are allowed. The most common
		// case we're worried about though is pentavalent
		// carbon. So we'll currently just throw an error
		// on that one case.
		
		if("C".equals(this.getSymbol()) && v>4){
			return true;
		}
		if("N".equals(this.getSymbol())){
			if(v==5){
				if(this.getCharge()==1){
					return false;
				}else{
					return true;
				}
			}
		}
		
		return false;
	}

	@Override
	public boolean isRGroupAtom() {
		return getRGroupIndex().isPresent();
	}

	@Override
	public boolean isPseudoAtom() {
		return CdkUtil.isPseudoAtom(atom);
	}

	@Override
	public OptionalInt getRGroupIndex() {
		
		return getPseudoAtomField(iPseudoAtom -> {
			String label = iPseudoAtom.getLabel();
			 if(label !=null && RGroupQuery.isValidRgroupQueryLabel(label)) {
				 return OptionalInt.of(Integer.parseInt(label.substring(1)));
			 }
			 return OptionalInt.empty();
		},
				OptionalInt::empty);

	}


	@Override
	public void setRGroup(Integer rGroup) {
		if(rGroup ==null || rGroup < 1) {
			if(isRGroupAtom()) {
				setAlias(null);
			}
		}else {
			setAlias("R"+rGroup);
		}
	}

	private <T> T getPseudoAtomField(Function<IPseudoAtom, T> function, Supplier<T> emptySupplier) {
		IAtom deref = AtomRef.deref(atom);
		if(deref instanceof IPseudoAtom) {
			IPseudoAtom iPseudoAtom = (IPseudoAtom)deref;
			return function.apply(iPseudoAtom);
		}
		return emptySupplier.get();
	}

	@Override
	public Optional<String> getAlias() {
		return getPseudoAtomField( a-> Optional.ofNullable(a.getLabel()), Optional::empty);
	}

	@Override
	public void setAlias(String alias) {
		changeToPseudoAtomIfNeeded((pseudoAtom, newObj) -> {
			pseudoAtom.setLabel(alias);
			if(newObj) {
				pseudoAtom.setSymbol(alias);
			}
			
		});
		
	}
	
	

	//borrowed from CDK MDLV20000Reader with minor adjustments
	//to make it take consumers for various operations
	
	private void changeToPseudoAtomIfNeeded(BiConsumer<IPseudoAtom, Boolean> consumer) {
		IAtomContainer container = parent.getContainer();
	 final IPseudoAtom pseudoAtom = atom instanceof IPseudoAtom ? (IPseudoAtom) atom : container.getBuilder()
             .newInstance(IPseudoAtom.class);
     if (atom.equals(pseudoAtom)) {
    	 consumer.accept(pseudoAtom, false);
     } else {
    	 //some consumers might change the symbol but by default reuse atom symbol
         pseudoAtom.setSymbol(atom.getSymbol());
         pseudoAtom.setAtomicNumber(0);
         pseudoAtom.setPoint2d(atom.getPoint2d());
         pseudoAtom.setPoint3d(atom.getPoint3d());
         pseudoAtom.setMassNumber(atom.getMassNumber());
         pseudoAtom.setFormalCharge(atom.getFormalCharge());
         pseudoAtom.setValency(atom.getValency());
         
         consumer.accept(pseudoAtom, true);
         //katzelda -chemkit note : need to update atom and bond caches
 		//Looking into the source code for this method
         //it doesn't make  new IBond objects just
         //updates the fields for  IBond objects
         //so we should  be OK and don't have to update our caches...
         
 		 // XXX: would be faster to track all replacements and do it all in one
         AtomContainerManipulator.replaceAtomByAtom(container, atom, pseudoAtom);
        
         this.atom = pseudoAtom;
     }
	}

	@Override
	public String toString() {
		return "CdkAtom [atom=" + atom + "]";
	}



	@Override
	public int getAtomIndexInParent() {
		return parent.getContainer().indexOf(atom);
	}

	
	

}
