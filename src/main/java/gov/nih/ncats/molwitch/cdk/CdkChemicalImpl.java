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

import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.vecmath.Tuple2d;

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
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
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

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.AtomCoordinates;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Bond.BondType;
import gov.nih.ncats.molwitch.SGroup.SGroupBracket;
import gov.nih.ncats.molwitch.SGroup.SGroupType;
import gov.nih.ncats.molwitch.Stereocenter;
import gov.nih.ncats.molwitch.BondTable;
import gov.nih.ncats.molwitch.ChemicalSource;
import gov.nih.ncats.molwitch.ChemkitException;
import gov.nih.ncats.molwitch.Chirality;
import gov.nih.ncats.molwitch.DoubleBondStereochemistry;
import gov.nih.ncats.molwitch.ExtendedTetrahedralChirality;
import gov.nih.ncats.molwitch.GraphInvariant;
import gov.nih.ncats.molwitch.SGroup;
import gov.nih.ncats.molwitch.TetrahedralChirality;
import gov.nih.ncats.molwitch.isotopes.Isotope;
import gov.nih.ncats.molwitch.spi.ChemicalImpl;
import gov.nih.ncats.common.util.CachedSupplier;

public class CdkChemicalImpl implements ChemicalImpl<CdkChemicalImpl>{


	private final ConcurrentHashMap<IAtom, CdkAtom> atoms = new ConcurrentHashMap<>();
	private final ConcurrentHashMap<IBond, CdkBond> bonds = new ConcurrentHashMap<>();
	
	private IAtomContainer container;
	
	private boolean isAromatic;
	
	CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
	CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(SilentChemObjectBuilder.getInstance());

    Map<Integer, Map<Integer, CdkBond>> bondMap = new HashMap<>();

    
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
	    	return null;
    }
    	);
    CachedSupplier<Void> ringsSearcherSupplier = CachedSupplier.of(()->{
		    	AllRingsFinder arf = new AllRingsFinder();
		    	IRingSet ringSet;
				try {
					ringSet = arf.findAllRings(container);
					for (int ir = 0; ir < ringSet.getAtomContainerCount(); ir++) {
			            IRing ring = (IRing) ringSet.getAtomContainer(ir);
			            for (int jr = 0; jr < ring.getAtomCount(); jr++) {
			                IAtom aring = ring.getAtom(jr);
			                aring.setFlag(CDKConstants.ISINRING, true);
			            }
			            
			            for (int jr = 0; jr < ring.getBondCount(); jr++) {
			                IBond aring = ring.getBond(jr);
			                aring.setFlag(CDKConstants.ISINRING, true);
			            }
					}
					return null;
				} catch (CDKException e) {
					throw new RuntimeException(e);
				}
		    	
		    	
		
    });
    
	private final ChemicalSource source;
	public CdkChemicalImpl(IAtomContainer container, Supplier<? extends ChemicalSource> source) {
		this(container, source.get());
	}
	public CdkChemicalImpl(IAtomContainer container, ChemicalSource source) {
		this.container = container;
		this.source = source;
        
		   for (IAtom atom : container.atoms()) {
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
		   }
		   
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
		IAtom atom = new org.openscience.cdk.Atom(symbol);
		   
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
		ringsSearcherSupplier.resetCache();
		cahnIngoldPrelogSupplier.resetCache();
		perceiveAtomTypesOfNonQueryAtoms.resetCache();
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
		
		return AtomContainerManipulator.getNaturalExactMass(container);
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

	
	@Override
	public void aromatize() {
		try {
			
//			fix();
			 new Aromaticity(ElectronDonation.daylight(),
                     Cycles.all())
			 .apply(container);
			isAromatic = true;
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private void percieveAtomTypeAndConfigureNonQueryAtoms() throws CDKException{
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(container.getBuilder());
        for (IAtom atom : container.atoms()) {
        		
            IAtomType matched = matcher.findMatchingAtomType(container, atom);
            if (matched != null) {
            		try{
            			AtomTypeManipulator.configure(atom, matched);
            		}catch(NullPointerException e) {
            			//ignore
            		}
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
			coordinateGenerator.generateCoordinates();
	
			container = coordinateGenerator.getMolecule();
		}catch(CDKException e) {
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
            else{
                kekulize();
            }
			if(options.computeCoords){
				StructureDiagramGenerator coordinateGenerator = new StructureDiagramGenerator(container);
				coordinateGenerator.generateCoordinates();

				container = coordinateGenerator.getMolecule();
			}

			if(options.computeStereo){
	            boolean has2dCoords=true;
	            boolean has3dCoords = true;
	            for(IAtom atom : container.atoms()){
	                if(atom.getPoint3d() ==null){
	                    has3dCoords=false;
	                }
	                if(atom.getPoint2d() ==null){
	                    has2dCoords=false;
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
		
		isAromatic=false;
//		fix();
		//setImplicitHydrogens();
	    try {
			Kekulization.kekulize(container);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
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
		for(IStereoElement se : container.stereoElements()){
			if(se instanceof ExtendedTetrahedral){
				
				 list.add(new CdkExtendedTetrahedralChirality((ExtendedTetrahedral)se));

			}
		}
		return list;
	}
	@Override
	public List<TetrahedralChirality> getTetrahedrals() {
		List<TetrahedralChirality> list = new ArrayList<>();

		cahnIngoldPrelogSupplier.get();
		
		for(IStereoElement se : container.stereoElements()){
			if(se instanceof ITetrahedralChirality){
				 ITetrahedralChirality chirality = (ITetrahedralChirality) se;
				 list.add(new CdkTetrahedralChirality(chirality));

			}
		}
		return list;
		
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
		//TODO CDK may deprecate this and add normal accessor in future version.
		List<Sgroup> sgroups = container.getProperty(CDKConstants.CTAB_SGROUPS);
		
		if(sgroups ==null) {
			return Collections.emptyList();
		}
		List<SGroup> chemkitGroups = new ArrayList<>();
		for(Sgroup group : sgroups) {
			chemkitGroups.add(new CDKSgroupAdapter(group));
		}
		return chemkitGroups;
	}
	
	
	@Override
	public void removeSGroup(SGroup sgroup) {
		List<Sgroup> sgroups = container.getProperty(CDKConstants.CTAB_SGROUPS);
		if(sgroups ==null) {
			return;
		}
		sgroups.remove( ((CDKSgroupAdapter)sgroup).sgroup);
		
	}
	@Override
	public SGroup addSgroup(SGroupType type) {
		SGroupType typeToUse = type==null? SGroupType.GENERIC : type;
		List<Sgroup> sgroups = container.getProperty(CDKConstants.CTAB_SGROUPS);
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
		Set<Sgroup> sgroups = container.getProperty(CDKConstants.CTAB_SGROUPS);
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
	        bond.getBegin().setFlag(CDKConstants.ISAROMATIC, true);
	        bond.getEnd().setFlag(CDKConstants.ISAROMATIC, true);
	        bond.setFlag(CDKConstants.ISAROMATIC, true);
		}
		container.addBond( bond);
		setImplicitHydrogens(((CdkAtom)atom1).getAtom());
		setImplicitHydrogens(((CdkAtom)atom2).getAtom());

		setDirty();
		return getCdkBondFor(bond);
	}
	protected void setImplicitHydrogens(){
		for(int i = 0; i< container.getAtomCount(); i++){
			setImplicitHydrogens(container.getAtom(i));
		}
	}
	protected int setImplicitHydrogens(IAtom atom){
		//compute it
		try {
			CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(SilentChemObjectBuilder.getInstance());
			 IAtomType type = matcher.findMatchingAtomType(container, atom);
			
			 AtomTypeManipulator.configure(atom, type);
			//includes explicit and implicit Hs
			 Integer neighborCount = atom.getFormalNeighbourCount();
			 if(neighborCount ==null){
				 neighborCount = 0;
			 }
			 int implicitH = neighborCount; 
			 for(IBond bond : container.getConnectedBondsList(atom)){
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
	private static Chirality convetCIPValueToChirality(IAtom center) {
		String value =center.getProperty(CDKConstants.CIP_DESCRIPTOR);
		if("S".equals(value)) {
			return Chirality.S;
		}
		if("R".equals(value)) {
			return Chirality.R;
		}
		if("NONE".equals(value)) {
			return Chirality.Non_Chiral;
		}
		System.out.println("stereo is not S or R but it's '" + value +"'");
		return Chirality.Unknown;
	}
	private static Chirality convertCdkStereoToChirality(Stereo stereoType) {
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
		public Atom getCenterAtom() {
			
			return getCdkAtomFor(chirality.getChiralAtom());
		}

		@Override
		public Chirality getChirality() {
			 return convetCIPValueToChirality(chirality.getChiralAtom());
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
		public boolean bracketsSupported() {
			return true;
		}

		@Override
		public Optional<String> getSubscript() {
			return Optional.ofNullable(sgroup.getSubscript());
		}

		@Override
		public Optional<String> getSuperscript() {
			return Optional.ofNullable(sgroup.getSubscript());
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
		return ret.toString();
	}

	


	@Override
	public void setProperty(String key, String value) {
		container.setProperty(key, value);
		
	}



	@Override
	public Iterator<Entry<String, String>> properties() {
		return new Iterator<Entry<String, String>>(){
			Iterator<Entry<Object,Object>> iter = container.getProperties().entrySet().iterator();

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public Entry<String, String> next() {
				final Entry<Object,Object> next = iter.next();
				
				return new Entry<String,String>() {
					
					@Override
					public String setValue(String value) {
						throw new UnsupportedOperationException();
					}
					
					@Override
					public String getValue() {
						Object v = next.getValue();
						if(v ==null){
							return null;
						}
						return v.toString();
					}
					
					@Override
					public String getKey() {
						return next.getKey().toString();
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
