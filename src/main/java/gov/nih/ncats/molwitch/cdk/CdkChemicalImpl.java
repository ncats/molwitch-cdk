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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.vecmath.Tuple2d;

import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.cip.CIPTool;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.sgroup.Sgroup;
import org.openscience.cdk.sgroup.SgroupBracket;
import org.openscience.cdk.sgroup.SgroupKey;
import org.openscience.cdk.sgroup.SgroupType;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.StereoElementFactory;
import org.openscience.cdk.stereo.Stereocenters;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;
import org.xmlcml.cml.base.CMLElement.Hybridization;

import gov.nih.ncats.common.util.CachedSupplier;
import gov.nih.ncats.common.util.CachedSupplierGroup;
import gov.nih.ncats.common.util.Unchecked;
import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.AtomCoordinates;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Bond.BondType;
import gov.nih.ncats.molwitch.BondTable;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.ChemicalSource;
import gov.nih.ncats.molwitch.ChemkitException;
import gov.nih.ncats.molwitch.Chirality;
import gov.nih.ncats.molwitch.DoubleBondStereochemistry;
import gov.nih.ncats.molwitch.ExtendedTetrahedralChirality;
import gov.nih.ncats.molwitch.GraphInvariant;
import gov.nih.ncats.molwitch.SGroup;
import gov.nih.ncats.molwitch.SGroup.SGroupBracket;
import gov.nih.ncats.molwitch.SGroup.SGroupType;
import gov.nih.ncats.molwitch.Stereocenter;
import gov.nih.ncats.molwitch.TetrahedralChirality;
import gov.nih.ncats.molwitch.isotopes.Isotope;
import gov.nih.ncats.molwitch.spi.ChemicalImpl;

public class CdkChemicalImpl implements ChemicalImpl<CdkChemicalImpl>{


	private final ConcurrentHashMap<IAtom, CdkAtom> atoms = new ConcurrentHashMap<>();
	private final ConcurrentHashMap<IBond, CdkBond> bonds = new ConcurrentHashMap<>();
	
	private IAtomContainer container;
	
	private boolean isAromatic;
	private static IChemObjectBuilder CHEM_OBJECT_BUILDER = SilentChemObjectBuilder.getInstance();
	private CDKHydrogenAdder hydrogenAdder;
//	CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(SilentChemObjectBuilder.getInstance());

    Map<Integer, Map<Integer, CdkBond>> bondMap = new HashMap<>();

    private final Aromaticity      aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));

    CachedSupplier<Void> perceiveAtomTypesOfNonQueryAtoms = CachedSupplier.of(()->{
    	try {
    		percieveAtomTypeAndConfigureNonQueryAtoms();
    	}catch(CDKException e) {
    		
    	}
    	return null;
    }
    );
    
    CachedSupplier<Void> cahnIngoldPrelogSupplier = CachedSupplier.of(()->{
	    	CIPTool.label(container);
	    	Stereocenters  centers   = Stereocenters.of(container);
	    	for (int i = 0; i < container.getAtomCount(); i++) {
	    		switch(centers.stereocenterType(i)){
					case Non:
						break;
					case Para:
					case Potential:
						break;
					case True:
						IAtom ai=container.getAtom(i);
						if(ai.getProperty(CDKConstants.CIP_DESCRIPTOR) == null){
							// we only really set it if it's tetrahedral CARBON
							// right now,
							boolean bailout=false;
							for(IBond ib:ai.bonds()){
								if(!Order.SINGLE.equals(ib.getOrder())){
									bailout=true;
									break;
								}								
							}
							
							if(bailout){
//								System.out.println("It's real but ...");
//								System.out.println(centers.elementType(i));
								
								continue;
							}
							if("N".equals(ai.getSymbol())){
								if(ai.getFormalCharge()==0){
									//only charged nitrogens allowed
									continue;
								}
							}
							
							
							
							container.getAtom(i).setProperty(CDKConstants.CIP_DESCRIPTOR, "EITHER");
							
						}
						break;
					default:
						break;
	    		}
	    	}
	    	return null;
    }
    	);
    CachedSupplier<Void> ringsSearcherSupplier = CachedSupplier.of(()->{
		    	AllRingsFinder arf = new AllRingsFinder();
		    	IRingSet ringSet;
				try {
					//max#Rings must not be > atomCount
					ringSet = arf.findAllRingsInIsolatedRingSystem(container, Math.min(container.getAtomCount(),12));
					for (int ir = 0; ir < ringSet.getAtomContainerCount(); ir++) {
			            IRing ring = (IRing) ringSet.getAtomContainer(ir);
			            for (int jr = 0; jr < ring.getAtomCount(); jr++) {
			                IAtom aring = ring.getAtom(jr);
			                aring.setIsInRing(true);
			            }
			            
			            for (int jr = 0; jr < ring.getBondCount(); jr++) {
			                IBond aring = ring.getBond(jr);
			                aring.setIsInRing(true);
			            }
					}
					return null;
				} catch (CDKException e) {
					throw new RuntimeException(e);
				}
		    	
		    	
		
    });
	private final CachedSupplierGroup cachedSupplierGroup = new CachedSupplierGroup();
	private final ChemicalSource source;
	public CdkChemicalImpl(IAtomContainer container, Supplier<? extends ChemicalSource> source) {
		this(container, source.get());
	}
	public CdkChemicalImpl(IAtomContainer container, ChemicalSource source) {

	    if(!(container instanceof IQueryAtomContainer)){
	        boolean isQuery =false;
	        for(IAtom a : container.atoms()){
	            if(AtomRef.deref(a) instanceof IQueryAtom){
	                isQuery = true;
	                break;
                }
            }
	        if(!isQuery){
                for(IBond a : container.bonds()){
                    if(BondRef.deref(a) instanceof IQueryBond){
                        isQuery = true;
                        break;
                    }
                }
            }
	        if(isQuery){
	            container = new QueryAtomContainer(container, CHEM_OBJECT_BUILDER);
            }
        }
		this.container = container;
		this.source = source;
		hydrogenAdder = CDKHydrogenAdder.getInstance(container.getBuilder());
		cachedSupplierGroup.add(ringsSearcherSupplier);
		cachedSupplierGroup.add(cahnIngoldPrelogSupplier);
		cachedSupplierGroup.add(perceiveAtomTypesOfNonQueryAtoms);

		   /*for (IAtom atom : container.atoms()) {
			   //query atoms don't have everything set
			   //so ignore them for now
			   if(atom instanceof IQueryAtom) {
				   continue;
			   }
		     IAtomType type;
			try {
				type = matcher.findMatchingAtomType(container, atom);
			} catch (CDKException e) {
				throw new IllegalStateException("error matching type", e);
			}

		     AtomTypeManipulator.configure(atom, type);
		   }*/
		   
	}
	

//	public void findRingAtoms(){
//		Atom
//	}
	
	@Override
	public void flipChirality(Stereocenter s) {
		for(Atom a : s.getPeripheralAtoms()) {
			for(Bond b : a.getBonds()) {
				gov.nih.ncats.molwitch.Bond.Stereo oldStereo = b.getStereo();
				gov.nih.ncats.molwitch.Bond.Stereo newStereo = oldStereo.flip();
				
				if(oldStereo !=newStereo) {
					IAtom center = CdkAtom.getIAtomFor(s.getCenterAtom());
					
				}
			}
		}
		
	}
	@Override
	public ChemicalSource getSource() {
		return source;
	}
	
	
	
	@Override
	//TODO: This is a weird thing to ask for ...
	// someone might want SSSR, but not typically the
	// smallest absolute ring.
	//
	// They may, however, want the smallest ring for
	// a particular bond/atom. Perhaps that's where
	// this came from?
	public int getSmallestRingSize() {
		IRingSet ringSet = Cycles.sssr(container).toRingSet();
		int numRings = ringSet.getAtomContainerCount();
		if(numRings ==0) {
			return 0;
		}
		int minSize=Integer.MAX_VALUE;
		for(int i=0; i< numRings; i++) {
			int s = ringSet.getAtomContainer(i).getAtomCount();
			if(s < minSize) {
				minSize= s;
			}
		}
		return minSize;
	}
	/**
	 * Get the {@link IAtomContainer}.
	 */
	@Override
	public Object getWrappedObject() {
		return container;
	}
	@Override
	public CdkChemicalImpl shallowCopy() {
		//shallow copy shares original atoms and bond objects
		return new CdkChemicalImpl(CdkUtil.getChemObjectBuilder().newInstance(IAtomContainer.class, container), source);
		
	}
	@Override
	public CdkChemicalImpl deepCopy() {
		
		try {
			return new CdkChemicalImpl(container.clone(), source);
		} catch (CloneNotSupportedException e) {
			throw new IllegalStateException();
		}
		
	}


	@Override
	public Atom addAtom(String symbol) {
		IAtom atom = CHEM_OBJECT_BUILDER.newAtom();
		atom.setSymbol(symbol);
		   
		return addAtom(atom);
	}

	


	@Override
	public void addChemical(ChemicalImpl<CdkChemicalImpl> other) {
		container.add((IAtomContainer) other.getWrappedObject());
		setDirty();
	}
	@Override
	public Atom addAtom(Atom a) {
		IAtom iatom = CdkAtom.getIAtomFor(a);
		return addAtom(iatom);
	}
	private Atom addAtom(IAtom iatom) {
		container.addAtom(iatom);
		setDirty();
		return getCdkAtomFor(iatom);
	}
	@Override
	public Atom addAtom(Isotope isotope) {
		IAtom a = CdkUtil.getChemObjectBuilder().newAtom();
		a.setSymbol(isotope.getSymbol());
		a.setAtomicNumber(isotope.getAtomicNumber());
		a.setExactMass(Double.valueOf(isotope.getMassNumber()));
		//TODO how to set the values with ranges/intervals?
		a.setNaturalAbundance(isotope.getRelativeAtomicMass().getValue().doubleValue());
		return addAtom(a);
	}
	
	@Override
	public Atom addAtomByAtomicNum(int atomicNumber) {
		return addAtom(PeriodicTable.getSymbol(atomicNumber));
		
	}

	private void setDirty() {

		cachedSupplierGroup.resetCache();
	}

	@Override
	public int getSGroupCount() {
		Collection<?> sgroups = getCdkSgroups();
		if(sgroups==null){
			return 0;
		}
		return sgroups.size();
	}

	@Override
	public void makeHydrogensExplicit() {
		
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(container);
		setDirty();
	}

	@Override
	public boolean hasImplicitHydrogens() {
		//short circuit as soon as you find one
		 int atomCount = container.getAtomCount();
	        for (int i = 0; i < atomCount; i++) {
	            Atom atom = getAtom(i);
	            Integer implicitNum = atom.getImplicitHCount();
	            if(implicitNum !=null && implicitNum.intValue() >0){
	            	return true;
	            }
	        }
	        return false;
	}


	@Override
	public void makeHydrogensImplicit() {
		AtomContainerManipulator.suppressHydrogens(container);
		setDirty();
	}
	public boolean isAromatic() {
		return isAromatic;
	}



	@Override
	public boolean hasCoordinates() {
		return getDim()>0;
	}
	@Override
	public boolean has2DCoordinates() {
		return getDim()==2;
	}
	@Override
	public boolean has3DCoordinates() {
		return getDim()==3;
	}
	
	private int getDim() {
		boolean has2dCoords=true;
        boolean has3dCoords = true;
        for(IAtom atom : container.atoms()){
            if(atom.getPoint3d() ==null){
                has3dCoords=false;
            }
            if(atom.getPoint2d() ==null){
                has2dCoords=false;
            }
            if(!has2dCoords && !has3dCoords) {
            	return 0;
            }
        }
        if(has3dCoords) {
        		return 3;
        }
        if(has2dCoords) {
    		return 2;
        }
        return 0;
	}
	@Override
	public String getName() {
		return container.getID();
	}

	@Override
	public void setName(String name) {
		container.setID(name);
		//some cdk objects like mol writers look for the title property
		container.setProperty(CDKConstants.TITLE, name);
	}

	CdkAtom getCdkAtomFor(IAtom atom){
        return atoms.computeIfAbsent(atom, a -> new CdkAtom(a, this));

	}
	
	CdkBond getCdkBondFor(IBond bond){

        //IBond doesn't implement equals or hashcode()!!!!
        //
        for(Entry<IBond, CdkBond> entry : bonds.entrySet()){
            IBond key = entry.getKey();
            if(  (bond.getAtom(0).equals(key.getAtom(0)) && bond.getAtom(1).equals(key.getAtom(1)) )
                || (bond.getAtom(1).equals(key.getAtom(0)) && bond.getAtom(0).equals(key.getAtom(1)))
            ){
                return entry.getValue();
            }
        }
		return bonds.computeIfAbsent(bond, b ->{
        //    System.out.println("creating new cdk wrapper bond for " + b);
            return  new CdkBond(b, this);
        });
	}
	
	List<CdkBond> getBondsFor(IAtom atom){
		List<CdkBond> bonds = new ArrayList<>();
		for(IBond bond : container.bonds()){
			if(bond.contains(atom)){
				bonds.add(getCdkBondFor(bond));
			}
		}
		return bonds;
	}

	@Override
	public double getMass() {
		perceiveAtomTypesOfNonQueryAtoms.get();
		boolean recomputeHydrogens=false;
		for(IAtom a : container.atoms()){
			if(a.getImplicitHydrogenCount()==null){
				recomputeHydrogens=true;
				break;
			}
//		setImplicitHydrogens()
		};
		if(recomputeHydrogens){
			setImplicitHydrogens();
		}
		return AtomContainerManipulator.getMass(container, AtomContainerManipulator.MolWeight);
	}

	@Override
	public int getAtomCount() {
		return container.getAtomCount();
	}

	@Override
	public int getBondCount() {
		return container.getBondCount();
	}

	@Override
	public Atom getAtom(int i) {
		return getCdkAtomFor(container.getAtom(i));
	}
	
	@Override
	public Atom removeAtom(int i){
		IAtom atom = container.getAtom(i);
		for(SGroup sgroup : getSGroups()){
			((CDKSgroupAdapter)sgroup).removeAtom(atom);
		}
		container.removeAtom(atom);
		Atom retAtom = getCdkAtomFor(atom);
		setDirty();
		return retAtom;
	}

	@Override
	public Atom removeAtom(Atom a){
		if(!(a instanceof CdkAtom)){
			throw new IllegalArgumentException("wrong type");
		}
		IAtom iatom = ((CdkAtom)a).getAtom();
		for(SGroup sgroup : getSGroups()){
			sgroup.removeAtom(a);

		}
		
		container.removeAtom(iatom);
		
		setDirty();
		return a;
	}
	
	

	@Override
	public Bond removeBond(int i) {
		IBond ibond = container.removeBond(i);
		Bond b = getCdkBondFor(ibond);
		setDirty();
		return b;
	}
	@Override
	public Bond removeBond(Bond b) {
		
		IBond iBond= ((CdkBond)b).getBond();
		container.removeBond(iBond);
		setDirty();
		return b;
	}
	
	
	@Override
	public Bond removeBond(Atom a, Atom b) {
		IBond iBond = container.removeBond(CdkAtom.getIAtomFor(a), CdkAtom.getIAtomFor(b));
		setDirty();
		return getCdkBondFor(iBond);
	}
	public IAtomContainer getContainer() {
		return container;
	}
	
	private <E extends Throwable> void doWithQueryFixes(Unchecked.ThrowingRunnable<E> r, boolean mustSet) throws E{
		Map<IAtom, Integer> oldAI = new HashMap<>();
		Map<IAtom, Integer> oldIH = new HashMap<>();
		Map<IAtom, Integer> oldFC = new HashMap<>();
		Map<IBond, Order> oldOrder = new HashMap<>();
//		
		for(IAtom atom : container.atoms()){
	         if(atom.getAtomicNumber()==null){
	        	 oldAI.put(atom, null);
	        	 atom.setAtomicNumber(2);
	         }
	         if(atom.getFormalCharge()==null){
	        	 oldFC.put(atom, null);
	        	 atom.setFormalCharge(0);
	         }
	         if(atom.getImplicitHydrogenCount()==null){
	        	 oldIH.put(atom, null);
	        	 //this is really unfortunate, because it 
	        	 atom.setImplicitHydrogenCount(0);	        	 
	         }
	         
	    } 
		for(IBond bond : container.bonds()){
			if(bond.getOrder()==null || bond.getOrder().equals(Order.UNSET)){
				Order oorder=bond.getOrder();
				if(bond instanceof QueryBond){
					QueryBond qb = (QueryBond)bond;
					if(qb.getExpression().type().equals(Expr.Type.ALIPHATIC_ORDER)||qb.getExpression().type().equals(Expr.Type.ORDER)){
						
						bond.setOrder(Order.values()[qb.getExpression().value()-1]);
						oldOrder.put(bond, oorder);
					}else{
						bond.setOrder(Order.SINGLE);
						oldOrder.put(bond, oorder);
					}
				}else{
					if(mustSet){
						bond.setOrder(Order.SINGLE);
						oldOrder.put(bond, oorder);
					}
				}
			}
		}
		
		try{
			r.run();
		}finally{
	        for(IAtom at: oldAI.keySet()){
	        	at.setAtomicNumber(null);
	        }
	        for(IAtom at: oldIH.keySet()){
	        	at.setImplicitHydrogenCount(null);
	        }
	        for(IAtom at: oldFC.keySet()){
	        	at.setFormalCharge(null);
	        }
	        for(IBond ib: oldOrder.keySet()){
	        	ib.setOrder(oldOrder.get(ib));
	        }
		}
	}

	
	@Override
	public void aromatize() {
		//CDK errors out if there are no atoms
		//so skip this function
		if(container.isEmpty()){
			isAromatic = true;
			return;
		}
		try {
			// This isn't pretty, but stereochemistry query bonds
			// end up messing up the aromaticity model
			// so they are explicitly set to be single or aromatic bonds
			// with implied order of 1. In practice, this may end up
			// causing problems, but really only when the bond
			// turns out not to be aromatic by the model ...
			// this will need to be evaluated.
			
			for(IBond bond : container.bonds()){
	            if(bond instanceof IQueryBond){
	            	IQueryBond iq = (IQueryBond)bond;
					QueryBond qb= (QueryBond)BondRef.deref(iq);
					Expr exp=qb.getExpression();
					CdkUtil.navNodes(exp, (e,l)->{
						if(l.type().equals(Expr.Type.STEREOCHEMISTRY)){
							l.setPrimitive(Expr.Type.SINGLE_OR_AROMATIC);
							iq.setOrder(Order.SINGLE);
						}
					},0);
	            }
	        }
			
			perceiveAtomTypesOfNonQueryAtoms.get();
//			fix();

			kekulize();

			setImplicitHydrogens();
			doWithQueryFixes(()->aromaticity.apply(container),true);
			
			
			// Due to the way queries are handled, an attempt
			// to aromatize atoms/bonds DOES mark them as aromatic,
			// but it doesn't update any underlying expression used
			// for matching. We need to look for cases where
			// a query expression explicitly asks for aliphatic 
			// bonds and elements and change them to be asking for
			// just the element and just the bond
			
			for(IAtom atom: container.atoms()){
				if(atom.isAromatic() && atom instanceof IQueryAtom){
					IQueryAtom iq = (IQueryAtom)atom;
					QueryAtom qa= (QueryAtom)AtomRef.deref(iq);
					Expr exp=qa.getExpression();
					CdkUtil.navNodes(exp, (e,l)->{
						if(l.type().equals(Expr.Type.ALIPHATIC_ELEMENT)){
							l.setPrimitive(Expr.Type.ELEMENT, l.value());
						}
					},0);	
				}	
			}
			for(IBond bond : container.bonds()){
	            if(bond instanceof IQueryBond){
	            	IQueryBond iq = (IQueryBond)bond;
					QueryBond qb= (QueryBond)BondRef.deref(iq);
					Expr exp=qb.getExpression();
					if(bond.isAromatic()){
						CdkUtil.navNodes(exp, (e,l)->{
							if(l.type().equals(Expr.Type.ALIPHATIC_ORDER)){
								l.setPrimitive(Expr.Type.IS_AROMATIC);
							}
						},0);
					}
	            }
	        }
			
			isAromatic = true;
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//		isAromatic = true;
	}
	
	private void percieveAtomTypeAndConfigureNonQueryAtoms() throws CDKException{
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(container.getBuilder());
        for (IAtom atom : container.atoms()) {

        	try{
            IAtomType matched = matcher.findMatchingAtomType(container, atom);
            if (matched != null) {
            		try{
            			AtomTypeManipulator.configure(atom, matched);
            		}catch(NullPointerException e) {
            			//ignore
            		}
            }
        	}catch(NullPointerException e) {
        		//ignore
        	}
        }
	}
	
	@Override
	public void expandSGroups() {
		//TODO
		
	}
	@Override
	public void generateCoordinates() throws ChemkitException{
		try {
			StructureDiagramGenerator coordinateGenerator = new StructureDiagramGenerator(container);
			doWithQueryFixes(coordinateGenerator::generateCoordinates,false);
			container = coordinateGenerator.getMolecule();
		}catch(Exception e) {
			throw new ChemkitException(e.getMessage(), e);
		}
		
	}
	@Override
	public void prepareForBuild(PreparationOptions options) {
		
		
		try {
			perceiveAtomTypesOfNonQueryAtoms.get();
			
			if(options.makeHydrogensExplicit){
				makeHydrogensExplicit();
			}else{
				hydrogenAdder.addImplicitHydrogens(container);
			}

            if(options.aromatize){
                aromatize();
            }
			if(options.computeCoords){
				StructureDiagramGenerator coordinateGenerator = new StructureDiagramGenerator(container);
				coordinateGenerator.generateCoordinates();

				container = coordinateGenerator.getMolecule();
			}

			if(options.computeStereo){
	            boolean has3dCoords = true;
	            for(IAtom atom : container.atoms()){
	                if(atom.getPoint3d() ==null){
	                    has3dCoords=false;
	                    break;
	                }
	            }
	            StereoElementFactory stereoElementFactory;
	            if(has3dCoords){
	                stereoElementFactory = StereoElementFactory.using3DCoordinates(container);
	            }else{
	                stereoElementFactory = StereoElementFactory.using2DCoordinates(container);
	            }
	
				List<IStereoElement> stereo = stereoElementFactory.createAll();
	
				container.setStereoElements(stereo);

			}

		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}



	static int aromStatus(IAtomContainer mol) {
	    int res = 0;
	    for (IAtom atom : mol.atoms()) {
	        if (atom.getImplicitHydrogenCount() == null && atom.getFlag(CDKConstants.ISAROMATIC)) {
	            // N and P are ambiguous
	            if (atom.getAtomicNumber() == 7 || atom.getAtomicNumber() == 15)
	                res = 2; //ambiguous
	            else if (res < 2)
	                res = 1;
	        }
	    }
	    return res;
	}
	@Override
	public void kekulize() {
		try {
			doWithQueryFixes(()->{
				
				Kekulization.kekulize(container);

				//kekulize doesn't touch the aromatic bond flags
				//so we want to I guess?
				for(IBond bond : container.bonds()){
					bond.setIsAromatic(false);
				}
				isAromatic=false;
			},false);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
       
		
	}




	@Override
	public String getFormula() {
		return getFormula(true);
	}
	
	
	@Override
	public String getFormula(boolean includeImplicitHydrogen) {
		// must manually generate formula
		// code based on MolecularFormulaManipulator#getMolecularFormula methods
		// in CDK

		IMolecularFormula formula = container.getBuilder().newInstance(IMolecularFormula.class);

		int charge = 0;
		IAtom hAtom = null;
		for (IAtom iAtom : container.atoms()) {
			formula.addIsotope(iAtom);
			if (iAtom.getFormalCharge() != null){
				charge += iAtom.getFormalCharge();
			}

			if (includeImplicitHydrogen && iAtom.getImplicitHydrogenCount() != null
					&& (iAtom.getImplicitHydrogenCount() > 0)) {
				if (hAtom == null){
					hAtom = container.getBuilder().newInstance(IAtom.class, "H");
				}
				formula.addIsotope(hAtom, iAtom.getImplicitHydrogenCount());
			}
		}
		formula.setCharge(charge);

		return MolecularFormulaManipulator.getString(formula);

	}
	
	
	

	@Override
	public List<ExtendedTetrahedralChirality> getExtendedTetrahedrals() {
		List<ExtendedTetrahedralChirality> list = new ArrayList<>();
		cahnIngoldPrelogSupplier.get();
		//TODO: this only finds defined stereocenters, it doesn't
		//find stereo centers that could be defined later
		for(IStereoElement se : container.stereoElements()){
			if(se instanceof ExtendedTetrahedral){
				
				 list.add(new CdkExtendedTetrahedralChirality((ExtendedTetrahedral)se));

			}
		}
		return list;
	}
	@Override
	public List<TetrahedralChirality> getTetrahedrals() {

		cahnIngoldPrelogSupplier.get();
		
		
		Chemical c= new Chemical(this);
		return c.atoms()
		 .filter(ca->{
			 switch(ca.getChirality()){
			case Non_Chiral:
				break;
			case Parity_Either:
			case R:
			case S:
				return true;
			case Unknown:
				break;
			default:
				break;
			 
			 }
			 return false;
		 })
		 .map(ca->new CdkTetrahedralChirality2(ca))
		 .collect(Collectors.toList());		
	}

	@Override
	public Bond addBond(Bond b) {
		IBond ibond = CdkBond.getIBondFor(b);
		container.addBond(ibond);
		setDirty();
		return getCdkBondFor(ibond);
	}

	//TODO SGroups are stored like this by molfile
	// container.setProperty(CDKConstants.CTAB_SGROUPS, sgroupCpyList);
	//where sgroupCpyList is a Set<SGroup> 
	@Override
	public List<SGroup> getSGroups(){
		List<Sgroup> sgroups = getCdkSgroups();


		if(sgroups ==null) {
			return Collections.emptyList();
		}
		List<SGroup> chemkitGroups = new ArrayList<>();
		for(Sgroup group : sgroups) {
			chemkitGroups.add(new CDKSgroupAdapter(group));
		}
		return chemkitGroups;
	}

	private List<Sgroup> getCdkSgroups() {
		//TODO CDK may deprecate this and add normal accessor in future version.
		return container.getProperty(CDKConstants.CTAB_SGROUPS);
	}


	@Override
	public void removeSGroup(SGroup sgroup) {
		List<Sgroup> sgroups = getCdkSgroups();
		if(sgroups ==null) {
			return;
		}
		sgroups.remove( ((CDKSgroupAdapter)sgroup).sgroup);
		
	}
	@Override
	public SGroup addSgroup(SGroupType type) {
		SGroupType typeToUse = type==null? SGroupType.GENERIC : type;
		List<Sgroup> sgroups = getCdkSgroups();
		if(sgroups ==null) {
			sgroups = new ArrayList<>();
			container.setProperty(CDKConstants.CTAB_SGROUPS, sgroups);
		}
		Sgroup cdkSgroup = new Sgroup();
		cdkSgroup.setType(SgroupType.parseCtabKey(typeToUse.getTypeName()));
		sgroups.add(cdkSgroup);
		return new CDKSgroupAdapter(cdkSgroup);
	}
	@Override
	public boolean hasSGroups() {
		List<Sgroup> sgroups = getCdkSgroups();
		return !(sgroups == null || sgroups.isEmpty());
	}
	@Override
	public Bond addBond(Atom atom1, Atom atom2, BondType type) {
		Order order = null;
		switch(type){
			case SINGLE : order = Order.SINGLE; break;
			case DOUBLE : order = Order.DOUBLE; break;
			case TRIPLE : order = Order.TRIPLE; break;
			case QUADRUPLE : order = Order.QUADRUPLE; break;
			case AROMATIC : order = Order.SINGLE; break; 
		}
		int index1 = indexOf(atom1);
		int index2 = indexOf(atom2);
		IBond bond = container.getBuilder().newInstance(IBond.class, 
					container.getAtom(index1), container.getAtom(index2), order);
		if(type == BondType.AROMATIC) {
	        bond.getBegin().setIsAromatic(true);
	        bond.getEnd().setIsAromatic(true);
	        bond.setIsAromatic(true);
		}
		container.addBond( bond);
		setImplicitHydrogens(((CdkAtom)atom1).getAtom());
		setImplicitHydrogens(((CdkAtom)atom2).getAtom());

		setDirty();
		return getCdkBondFor(bond);
	}
	protected void setImplicitHydrogens(){
		try {
			hydrogenAdder.addImplicitHydrogens(container);
		} catch (CDKException e) {
			e.printStackTrace();
		}
	}
	protected int setImplicitHydrogens(IAtom atom){
		//compute it
		try {
			CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(CHEM_OBJECT_BUILDER);
			 IAtomType type = matcher.findMatchingAtomType(container, atom);
			
			 AtomTypeManipulator.configure(atom, type);
			//includes explicit and implicit Hs
			 Integer neighborCount = atom.getFormalNeighbourCount();
			 if(neighborCount ==null){
				 neighborCount = 0;
			 }
			 int implicitH = neighborCount;
//			 System.out.println(implicitH + "  " + atom);
			 for(IBond bond : container.getConnectedBondsList(atom)){
//			 	System.out.println("\t"+bond.getOrder());
				 switch(bond.getOrder()){
				 	case SINGLE : implicitH--; break;
				 	case DOUBLE : implicitH-=2; break;
				 	case TRIPLE : implicitH-=3; break;
				 	case QUADRUPLE :implicitH-=4; break;
				 	case QUINTUPLE :implicitH-=5; break;
				 	case SEXTUPLE :implicitH-=6; break;
				 	default : // ? TODO do we subtract 1?
				 }
			 }
//			 System.out.println(implicitH);
			atom.setImplicitHydrogenCount(Math.max(0, implicitH));
			return implicitH;
			 
		}catch(CDKException e){
			throw new RuntimeException("error computing implicit H count for atom " +atom, e);
		}
	}
	
	@Override
	public Iterator<CdkChemicalImpl> connectedComponents(){
		Iterator<IAtomContainer> iter = ConnectivityChecker.partitionIntoMolecules(container).atomContainers().iterator();
		return new Iterator<CdkChemicalImpl>(){

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public CdkChemicalImpl next() {
				return new CdkChemicalImpl(iter.next(), getSource());
			}
			
		};
	}
	
	
	@Override
	public int indexOf(Atom a) {
		return container.indexOf(CdkAtom.getIAtomFor(a));
	}

	@Override
	public int indexOf(Bond b) {
		return container.indexOf(CdkBond.getIBondFor(b));
	}


	@Override
	public GraphInvariant getGraphInvariant() {
		int[][] g = GraphUtil.toAdjList(container);
		long[] inv =Canon.basicInvariants(container, g);
		return new CdkGraphInvariant(inv);
//		try {
//			return new CdkGraphInvariant(InChINumbersTools.(container));
//		} catch (CDKException e) {
//			e.printStackTrace();
//			//TODO should we throw an exception?
//			return null;
//		}
	}


	@Override
	public Bond getBond(int i) {		
		return getCdkBondFor(container.getBond(i));
	}

	@Override
	public BondTable getBondTable() {
		EdgeToBondMap bondMap = EdgeToBondMap.withSpaceFor(container);
		//bondMap is empty! need to populate it
		
		//only way appears to be to use adjList
		//then throw out the adjacency list
		//this is all because bondMap.put(i, j, bond) is private!
		
		GraphUtil.toAdjList(container, bondMap);
		
		return new CdkBondTable(bondMap);
	}

	private class CdkBondTable implements BondTable{

		private final EdgeToBondMap bondMap ;
		
		public CdkBondTable(EdgeToBondMap bondMap) {
			this.bondMap = bondMap;
		}

		@Override
		public Bond getBond(int i, int j) {
			return CdkChemicalImpl.this.getCdkBondFor(bondMap.get(i, j));
		}

		@Override
		public boolean bondExists(int i, int j) {
			return bondMap.get(i, j) !=null;
		}

		@Override
		public int getAtomCount() {
			return CdkChemicalImpl.this.getAtomCount();
		}
		
		
	}

	interface CDKStereocenter{
		public Set<Stereocenters.Stereocenter> removeStereoCenterFrom(IAtomContainer container);
		public Stereocenter flip();
	}
	private class CdkExtendedTetrahedralChirality implements ExtendedTetrahedralChirality, CDKStereocenter{
		
		private List<Atom> terminalAtoms;
		private List<Atom>  peripheralAtoms;
		private Atom centerAtom;
		private Chirality stereoType;
		
		public CdkExtendedTetrahedralChirality(ExtendedTetrahedral extendedTetrahedral) {
			
			terminalAtoms = toCdkAtomArray(extendedTetrahedral.findTerminalAtoms(container));
			//if a peripheral atom is implicit then the terminal atom is listed instead
			//for inchi calculations we have to take that into account
			peripheralAtoms = toCdkAtomArray(extendedTetrahedral.peripherals());
			
			centerAtom = getCdkAtomFor(extendedTetrahedral.focus());
			
			stereoType = convertCdkStereoToChirality(extendedTetrahedral.winding());
		}

		@Override
		public boolean isDefined() {
			//check if any terminal atoms have wedge or hash?
			for(Atom a : terminalAtoms){

				for(Bond b : a.getBonds()){
					//make sure it's to the peripheral atom
					if(peripheralAtoms.contains(b.getOtherAtom(a)) && b.getStereo() !=Bond.Stereo.NONE){
						return true;
					}
				}
			}
			return false;
		}

		@Override
		public Set<Stereocenters.Stereocenter> removeStereoCenterFrom(
				IAtomContainer container) {
			// TODO Auto-generated method stub
			return null;
		}



		@Override
		public Stereocenter flip() {
			// TODO Auto-generated method stub
			return null;
		}



		private List<Atom> toCdkAtomArray(IAtom[] as) {
			 return Collections.unmodifiableList(Arrays.stream(as)
				.map(CdkChemicalImpl.this::getCdkAtomFor)
				.collect(Collectors.toList()));
		}
		@Override
		public Chirality getChirality() {
			return stereoType;
		}

		/* (non-Javadoc)
		 * @see gov.nih.ncats.chemkit.spi.cdk.ExtendedTetrahedralChirality#getTerminalAtoms()
		 */
		@Override
		public List<Atom> getTerminalAtoms() {
			return terminalAtoms;
		}

		/* (non-Javadoc)
		 * @see gov.nih.ncats.chemkit.spi.cdk.ExtendedTetrahedralChirality#getPeripheralAtoms()
		 */
		@Override
		public List<Atom> getPeripheralAtoms() {
			return peripheralAtoms;
		}

		/* (non-Javadoc)
		 * @see gov.nih.ncats.chemkit.spi.cdk.ExtendedTetrahedralChirality#getCenterAtom()
		 */
		@Override
		public Atom getCenterAtom() {
			return centerAtom;
		}
		
		
	}
	private static Chirality convertCIPValueToChirality(IAtom center) {
		String value =center.getProperty(CDKConstants.CIP_DESCRIPTOR);
		if("S".equals(value)) {
			return Chirality.S;
		}
		if("R".equals(value)) {
			return Chirality.R;
		}
		if("EITHER".equals(value)) {
			return Chirality.Parity_Either;
		}
		if("NONE".equals(value)) {
			return Chirality.Non_Chiral;
		}
//		System.out.println("stereo is not S or R but it's '" + value +"'");
		return Chirality.Unknown;
	}
	private static Chirality convertCdkStereoToChirality(Stereo stereoType) {
		//I don't think this is right
		switch(stereoType){
			 case ANTI_CLOCKWISE : return Chirality.S;
			 case CLOCKWISE : return Chirality.R;
		 }
		 return null;
	}
	
	private class CdkTetrahedralChirality implements TetrahedralChirality{
		private final ITetrahedralChirality chirality;
		private final IAtom[] ligands;
		public CdkTetrahedralChirality(ITetrahedralChirality chirality) {
			this.chirality = chirality;
			ligands = chirality.getLigands();
		}

		@Override
		public boolean isDefined() {
			//let's just assume if any hash or wedge bonds
			//are set then it's defined...
			for(IBond b : chirality.getChiralAtom().bonds()){
				if(b.getStereo() != IBond.Stereo.NONE){
					return true;
				}
			}
			return false;
		}

		@Override
		public Atom getCenterAtom() {
			
			return getCdkAtomFor(chirality.getChiralAtom());
		}

		@Override
		public Chirality getChirality() {
			 return convertCIPValueToChirality(chirality.getChiralAtom());
		}

		

		@Override
		public Atom getLigand(int i) {
			return getCdkAtomFor(ligands[i]);
		}

		@Override
		public List<Atom> getPeripheralAtoms() {
			List<Atom> atoms = new ArrayList<Atom>(4);
			atoms.add(getLigand(0));
			atoms.add(getLigand(1));
			atoms.add(getLigand(2));
			atoms.add(getLigand(3));
			return atoms;
		}
	}
	
	private class CdkTetrahedralChirality2 implements TetrahedralChirality{
		private final Atom ia;
		private final Atom[] ligands;
		
		public CdkTetrahedralChirality2(Atom ai) {
			this.ia = ai;
			ligands = new Atom[4];
			List<Atom> nats= getCenterAtom().getNeighbors();
			
			for(int i=0;i<nats.size();i++){
				ligands[i] = nats.get(i);
			}
		}

		@Override
		public boolean isDefined() {
			Atom ca= getCenterAtom();
			Chirality c=ca.getChirality();
			// Chirality.R.equals(c) || Chirality.S.equals(c)
			// TODO: really, it can be defined and not be R/S,
			// as in some rings. Should fix.
			
			return Chirality.R.equals(c) || Chirality.S.equals(c);
		}

		@Override
		public Atom getCenterAtom() {
			return ia;
		}

		@Override
		public Chirality getChirality() {
			 return getCenterAtom().getChirality();
		}		

		@Override
		public Atom getLigand(int i) {
			return ligands[i];
		}

		@Override
		public List<Atom> getPeripheralAtoms() {
			List<Atom> atoms = new ArrayList<Atom>(4);
			atoms.add(getLigand(0));
			atoms.add(getLigand(1));
			atoms.add(getLigand(2));
			atoms.add(getLigand(3));
			return atoms;
		}
	}
	
	private class CDKSgroupAdapter implements SGroup{


		private Sgroup sgroup;
		
		public CDKSgroupAdapter( Sgroup sgroup) {
			this.sgroup = sgroup;
		}

		@Override
		public SGroupType getType() {
			return SGroupType.valueByTypeName(sgroup.getType().getKey());
			
		}

		@Override
		public Optional<String> getSuperatomLabel() {
			return Optional.empty();
		}

		@Override
		public Stream<Atom> getAtoms() {
			return sgroup.getAtoms().stream().map(a-> getCdkAtomFor(a));
		}

		@Override
		public Stream<Bond> getBonds() {
			return sgroup.getBonds().stream().map(b-> getCdkBondFor(b));
		}
		
		

		@Override
		public boolean bracketsTrusted() {
			return true;
		}

		@Override
		public SGroupConnectivity getConnectivity() {
			String value = sgroup.getValue(SgroupKey.CtabConnectivity);
			return SGroupConnectivity.parse(value);
		}

		@Override
		public Stream<SGroup> getParentHierarchy() {
			
			//CDK sgroup should never have nulls
			return sgroup.getParents().stream().map(CDKSgroupAdapter::new);
		}

		@Override
		public Stream<Atom> getOutsideNeighbors() {
			Set<IAtom> sgroupAtoms = sgroup.getAtoms();
			Set<IAtom> outside = new HashSet<>();
			for(IAtom a : sgroupAtoms) {
				for(IAtom a2 : container.getConnectedAtomsList(a)) {
					if(!sgroupAtoms.contains(a2)) {
						outside.add(a2);
					}
				}
			}
			return outside.stream().map(a-> getCdkAtomFor(a));
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + ((sgroup == null) ? 0 : sgroup.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			if (!(obj instanceof CDKSgroupAdapter)) {
				return false;
			}
			CDKSgroupAdapter other = (CDKSgroupAdapter) obj;
			if (!getOuterType().equals(other.getOuterType())) {
				return false;
			}
			if (sgroup == null) {
				if (other.sgroup != null) {
					return false;
				}
			} else if (!sgroup.equals(other.sgroup)) {
				return false;
			}
			return true;
		}

		private CdkChemicalImpl getOuterType() {
			return CdkChemicalImpl.this;
		}

		@Override
		public void addAtom(Atom a) {
			sgroup.addAtom(CdkAtom.getIAtomFor(a));
			
		}

		@Override
		public void addBond(Bond b) {
			sgroup.addBond(CdkBond.getIBondFor(b));
			
		}

		@Override
		public void removeAtom(Atom a) {
			sgroup.removeAtom(CdkAtom.getIAtomFor(a));
			
		}
		public void removeAtom(IAtom a) {
			sgroup.removeAtom(a);

		}

		@Override
		public void removeBond(Bond b) {
			sgroup.removeBond(CdkBond.getIBondFor(b));
			
		}

		@Override
		public PolymerSubType getPolymerSubType() {
			return PolymerSubType.valueByCode(sgroup.getValue(SgroupKey.CtabSubType));
		}

		@Override
		public void setPolymerSubType(PolymerSubType polymerSubtype) {
			if(polymerSubtype ==null) {
				sgroup.getAttributeKeys().remove(SgroupKey.CtabSubType);
			}else {
				sgroup.putValue(SgroupKey.CtabSubType, polymerSubtype.getCode());
			}
			
		}

		@Override
		public boolean hasBrackets() {
			List<SgroupBracket> brackets = sgroup.getValue(SgroupKey.CtabBracket);
			return !(brackets ==null || brackets.isEmpty());
		}

		@Override
		public List<SGroupBracket> getBrackets() {
			List<SgroupBracket> brackets = sgroup.getValue(SgroupKey.CtabBracket);
			if(brackets ==null){
				return Collections.emptyList();
			}
			return brackets.stream().map(SGroupBracketImpl::new).collect(Collectors.toList());
		}

        @Override
        public Optional<String> getSruLabel() {
			//NCATS/FDA often use SMT as the SRU label which CDK reads as subscript property
		    return getSubscript();
        }

        @Override
		public boolean bracketsSupported() {
			return true;
		}

		@Override
		public Optional<String> getSubscript() {
		    /*
		    public static final int ST_SUPERATOM = 0;
    public static final int ST_MULTIPLE = 1;
    public static final int ST_SRU = 2;
    public static final int ST_MONOMER = 3;
    public static final int ST_MER = 4;
    public static final int ST_COPOLYMER = 5;
    public static final int ST_CROSSLINK = 6;
    public static final int ST_MODIFICATION = 7;
    public static final int ST_MIXTURE = 8;
    public static final int ST_FORMULATION = 9;
    public static final int ST_DATA = 10;
    public static final int ST_ANY = 11;
    public static final int ST_GENERIC = 12;
    public static final int ST_COMPONENT = 13;

    public static final int SST_ALTERNATING = 1;
    public static final int SST_RANDOM = 2;
    public static final int SST_BLOCK = 3;
		    public String getSubscript() {
        if (this.sgroupType == ST_COMPONENT) {
            return this.sgroupSubscript == null ? "c" : this.sgroupSubscript;
        } else if (this.sgroupType == ST_MIXTURE) {
            return "mix";
        } else if (this.sgroupType == ST_FORMULATION) {
            return "f";
        } else if (this.sgroupSubscript != null) {
            return this.sgroupSubscript;
        } else if (this.sgroupType == ST_DATA) {
            return "";
        } else if (this.sgroupSubType == SST_RANDOM) {
            return "ran";
        } else {
            return this.sgroupType == ST_SRU ? "n" : "";
        }
    }
		     */

		    //try to match old jchem
            String actualSubscriptValue = sgroup.getSubscript();
            SGroupType type = getType();
            if(type == SGroupType.MULTIPLE){
                //subscript is the multiplier
                return Optional.ofNullable(actualSubscriptValue);
            }
            if(type==SGroupType.COMPONENT){
                if(actualSubscriptValue==null){
                    return Optional.of("c");
                }
                return Optional.of(actualSubscriptValue);
            }
            if(type == SGroupType.MIXTURE){
                return Optional.of("mix");
            }
            if(type == SGroupType.FORMULATION){
                return Optional.of("f");
            }
            if(actualSubscriptValue !=null){
                return Optional.of(actualSubscriptValue);
            }
            if(type == SGroupType.DATA){
                return Optional.empty();
            }
            if(getPolymerSubType() == PolymerSubType.RANDOM){
                return Optional.of("ran");
            }
            if(type == SGroupType.SRU){
                return Optional.of("n");
            }
			return Optional.empty();
		}

		@Override
		public Optional<String> getSuperscript(){
            SGroupType type = getType();
            if(type == SGroupType.SRU || type == SGroupType.MULTIPLE){

                //ignore it
                return Optional.empty();
            }
         SGroupConnectivity connectivity= getConnectivity();

         if(connectivity==null){
             return Optional.empty();
         }
         switch(connectivity){
             case HEAD_TO_HEAD: return Optional.of("hh");
             case HEAD_TO_TAIL: return Optional.of("ht");
             default : return Optional.empty();
         }
		}
		
		

	}


	@Override
	public List<DoubleBondStereochemistry> getDoubleBondStereochemistry() {
		List<DoubleBondStereochemistry> list = new ArrayList<>();
		for(IStereoElement se : container.stereoElements()){
			if(se instanceof IDoubleBondStereochemistry){
				IDoubleBondStereochemistry dbStereo = (IDoubleBondStereochemistry) se;
				list.add(new DoubleBondSeteochemistryImpl(dbStereo));
			}
		}
		return list;
	}
	
	private class SGroupBracketImpl implements SGroupBracket{

		private final SgroupBracket bracket;
		
		public SGroupBracketImpl(SgroupBracket bracket) {
			this.bracket = bracket;
		}

		@Override
		public AtomCoordinates getPoint1() {
			return toCoordinates(bracket.getFirstPoint());
		}

		@Override
		public AtomCoordinates getPoint2() {
			return toCoordinates(bracket.getSecondPoint());
		}
		
		
	}
	private static AtomCoordinates toCoordinates(Tuple2d tuple2d){
		return AtomCoordinates.valueOf(tuple2d.x, tuple2d.y);
	}
	
	private class DoubleBondSeteochemistryImpl implements DoubleBondStereochemistry{
		private final IBond stereoBond, doubleBond;
		private final IAtom ligands[];
		private final DoubleBondStereo stereo;
		private final IDoubleBondStereochemistry dbStereo;
		
		public DoubleBondSeteochemistryImpl(IDoubleBondStereochemistry dbStereo) {
			this.dbStereo = dbStereo;
			
			IBond[] surroundingBonds = dbStereo.getBonds();
            if (surroundingBonds[0] == null || surroundingBonds[1] == null)
                throw new RuntimeException("Cannot generate an InChI with incomplete double bond info");
            IDoubleBondStereochemistry.Conformation stereoType = dbStereo.getStereo();
            
            if(stereoType == IDoubleBondStereochemistry.Conformation.TOGETHER){
            	stereo = DoubleBondStereo.Z_CIS;
            }else if (stereoType == IDoubleBondStereochemistry.Conformation.OPPOSITE){
            	stereo = DoubleBondStereo.E_TRANS;
            }else{
                stereo = DoubleBondStereo.E_OR_Z;
            }
            
            this.stereoBond = dbStereo.getStereoBond();

            ligands = new IAtom[4];
            
            if (stereoBond.contains(surroundingBonds[0].getAtom(0))) {
                // first atom is A
            	ligands[1] = surroundingBonds[0].getAtom(0);
                ligands[0] = surroundingBonds[0].getAtom(1);
            } else {
                // first atom is X
            	ligands[0] = surroundingBonds[0].getAtom(0);
            	ligands[1] = surroundingBonds[0].getAtom(1);
            }
            if (stereoBond.contains(surroundingBonds[1].getAtom(0))) {
                // first atom is B
            	ligands[2] = surroundingBonds[1].getAtom(0);
            	ligands[3] = surroundingBonds[1].getAtom(1);
            } else {
                // first atom is Y
            	ligands[2] = surroundingBonds[1].getAtom(1);
            	ligands[3] = surroundingBonds[1].getAtom(0);
            }

            /*
            Stereochemistry specification for double bond stereochemistry.
            The data model defines the double atoms and two ligands attached to those two atoms,
            linearly connected with the double bond in the middle. The IBonds that define the stereo element
             are defined as an array where the bonds are sorted according to the linear connectivity.
              Thus, the first and third bonds are the two bonds attached on either side of the double bond,
              and the second bond is the double bond. The stereo annotation then indicates if the
               ligand atoms are in the cis positio
             */
            doubleBond =  container.getBond(stereoBond.getAtom(0), stereoBond.getAtom(1));
		}
		
		

        @Override
		public Bond getLigandBond(int i) {
        		return getCdkBondFor(dbStereo.getBonds()[i]);
		}



		private IBond findDoubleBond(IDoubleBondStereochemistry dbStereo) {
            for(IBond bond : dbStereo.getBonds()){
                if(bond.getOrder() == Order.DOUBLE){
                    return bond;
                }
            }
            throw new IllegalStateException("could not find double bond!!!");
        }

        @Override
		public DoubleBondStereo getStereo() {
			return stereo;
		}

		@Override
		public Bond getDoubleBond() {
            /*
            Stereochemistry specification for double bond stereochemistry.
            The data model defines the double atoms and two ligands attached to those two atoms,
            linearly connected with the double bond in the middle. The IBonds that define the stereo element
             are defined as an array where the bonds are sorted according to the linear connectivity.
              Thus, the first and third bonds are the two bonds attached on either side of the double bond,
              and the second bond is the double bond. The stereo annotation then indicates if the
               ligand atoms are in the cis positio
             */
            return getCdkBondFor(doubleBond);
		}

		@Override
		public Atom getLigand(int i) {
			return getCdkAtomFor(ligands[i]);
		}
		
	}


	@Override
	public String getProperty(String key) {
		Object ret= container.getProperty(key);
		if(ret==null){
			return null;
		}
		if(key.equals(CDKConstants.CTAB_SGROUPS)){
			System.out.println("here!!!!!!");
			new Exception("getting sgroup as string").printStackTrace();
		}
		return ret.toString();
	}


	@Override
	public void removeProperty(String name) {
		container.removeProperty(name);
	}

	@Override
	public void setProperty(String key, String value) {
		container.setProperty(key, value);
		
	}



	@Override
	public Iterator<Entry<String, String>> properties() {
		return new Iterator<Entry<String, String>>(){
			Iterator<Entry<Object,Object>> iter = container.getProperties().entrySet().iterator();

			Entry<Object,Object> next = getNextRealProperty();

			@Override
			public boolean hasNext() {
				return next !=null;
			}

			private Entry<Object, Object> getNextRealProperty(){
				while(iter.hasNext()){
					Entry<Object, Object> next = iter.next();
					//CDK stores some objects like Sgroups as properties we don't want to expose those
					if(CDKConstants.CTAB_SGROUPS.equals(next.getKey())){
						//don't include this!!
					}else{
						return next;
					}
				}
				return null;
			}
			@Override
			public Entry<String, String> next() {
				if(!hasNext()){
					throw new NoSuchElementException();
				}

				final Entry<Object,Object> retNext = this.next;
				this.next=getNextRealProperty();
				return new Entry<String,String>() {
					
					@Override
					public String setValue(String value) {
						throw new UnsupportedOperationException();
					}
					
					@Override
					public String getValue() {
						Object v = retNext.getValue();
						if(v ==null){
							return null;
						}
						return v.toString();
					}
					
					@Override
					public String getKey() {
						return retNext.getKey().toString();
					}
				};
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
				
			}
			
			
		};
	}
	
	
	

}
