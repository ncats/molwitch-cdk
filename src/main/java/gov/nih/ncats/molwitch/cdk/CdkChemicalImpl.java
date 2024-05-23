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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
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
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.vecmath.Tuple2d;

import org.junit.Assert;
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.cip.CIPTool;
import org.openscience.cdk.geometry.cip.CIPToolMod;
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

import gov.nih.ncats.common.Tuple;
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
import gov.nih.ncats.molwitch.Chirality;
import gov.nih.ncats.molwitch.DoubleBondStereochemistry;
import gov.nih.ncats.molwitch.ExtendedTetrahedralChirality;
import gov.nih.ncats.molwitch.GraphInvariant;
import gov.nih.ncats.molwitch.MolwitchException;
import gov.nih.ncats.molwitch.SGroup;
import gov.nih.ncats.molwitch.SGroup.SGroupBracket;
import gov.nih.ncats.molwitch.SGroup.SGroupType;
import gov.nih.ncats.molwitch.Stereocenter;
import gov.nih.ncats.molwitch.TetrahedralChirality;
import gov.nih.ncats.molwitch.isotopes.Isotope;
import gov.nih.ncats.molwitch.isotopes.NISTIsotopeFactory;
import gov.nih.ncats.molwitch.spi.ChemicalImpl;

import static org.openscience.cdk.geometry.cip.CIPToolMod.getSizeOfLargestRingSystem;

public class CdkChemicalImpl implements ChemicalImpl<CdkChemicalImpl>{


	private final ConcurrentHashMap<IAtom, CdkAtom> atoms = new ConcurrentHashMap<>();
	private final ConcurrentHashMap<IBond, CdkBond> bonds = new ConcurrentHashMap<>();
	
	private IAtomContainer container;
	
	private boolean isAromatic;
	private static IChemObjectBuilder CHEM_OBJECT_BUILDER = SilentChemObjectBuilder.getInstance();
	private CDKHydrogenAdder hydrogenAdder;

	//this parameter allows us to control whether to use the newer (more accurate but prohibitively slow in some complex
	// chemical systems) or the older (faster but incorrect, in some cases) methods of atom labelling.
	// The default value is 7, a guess and a setter allows a calling routine to change this cutoff
	private static int complexityCutoff = 7;

	public static int getComplexityCutoff() {
		return CdkChemicalImpl.complexityCutoff;
	}

	public static void setComplexityCutoff(int complexityCutoff) {
		CdkChemicalImpl.complexityCutoff = complexityCutoff;
	}

	public static Integer getMaxUndefinedStereoCenters() {
		return CdkChemicalImpl.maxUndefinedStereoCenters;
	}

	public static void setMaxUndefinedStereoCenters(int maxUndefinedStereoCenters) {
		CdkChemicalImpl.maxUndefinedStereoCenters = maxUndefinedStereoCenters;
	}

	//This parameter allows us to prevent a long calculation that enumerates all possibilities of undefined stereocenters.
	//  this should typically be set to values < 10
	private static int maxUndefinedStereoCenters = 5;

//	CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(SilentChemObjectBuilder.getInstance());

    Map<Integer, Map<Integer, CdkBond>> bondMap = new HashMap<>();

    private final Aromaticity      aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.all(6)));
    
    private final static int[] mostStable = new int[119];
    private final static BitSet cdkMissing = new BitSet();

	private boolean deepChirality = true;

	static{
		
		mostStable[94]=244;
		mostStable[96]=247;
		mostStable[43]=97; 
		mostStable[93]=237;
		mostStable[91]=231;
		mostStable[95]=243;
		mostStable[88]=226;
		mostStable[97]=247;
		mostStable[98]=251;
		mostStable[84]=209;
		mostStable[89]=227;
		mostStable[61]=145;
		mostStable[99]=252;
		mostStable[100]=257;
		mostStable[101]=258;
		mostStable[86]=222;
		mostStable[105]=268;
		mostStable[103]=266;
		mostStable[85]=210;
		mostStable[104]=267;
		mostStable[102]=259;
		mostStable[87]=223;
		mostStable[106]=269;
		mostStable[111]=282;
		mostStable[107]=270;
		mostStable[112]=285;
		mostStable[108]=269;
		mostStable[110]=281;
		mostStable[113]=286;
		mostStable[109]=278;
		mostStable[114]=289;
		mostStable[115]=290;
		mostStable[116]=293;
		mostStable[117]=294;
		mostStable[118]=294;
		
		cdkMissing.set(43);
		cdkMissing.set(61);
		for(int i=84;i<=89;i++)cdkMissing.set(i);
		for(int i=93;i<=118;i++)cdkMissing.set(i);
		
    }

	CachedSupplier<Boolean> complexitySupplier =CachedSupplier.of(()->{
		int sizeOfLargestRingSystem = getSizeOfLargestRingSystem(this);
		Logger.getLogger(this.getClass().getName()).info("getSizeOfLargestRingSystem(this): " + sizeOfLargestRingSystem);
		return sizeOfLargestRingSystem > complexityCutoff;
	});

	CachedSupplier<Integer> perceiveAtomTypesOfNonQueryAtoms = CachedSupplier.of(()->{
    	try {
    		percieveAtomTypeAndConfigureNonQueryAtoms();
    	}catch(CDKException e) {
    		
    	}
    	
    	return 1;
    }
    );


	CachedSupplier<Integer> cahnIngoldPrelogSupplier = CachedSupplier.of(()->{
		boolean isComplex = complexitySupplier.get();
        try {
            makeStereoElms() ;
        }catch(Exception e) {
            //TODO: log exception here. It usually happens because of query atoms
        }
            withModifiedForm(c->{
               CdkChemicalImpl cimp = (CdkChemicalImpl) c.getImpl();
//               cimp.
               
               //Due to a bug in CDK, this doesn't work as expected
               //CIPTool.label(cimp.getContainer());
               //We currently have to use a modified form for now
				if( isComplex) {
					Logger.getLogger(this.getClass().getName()).fine("This molecule is considered complex");
					CIPTool.label(cimp.getContainer());
				} else {
					Logger.getLogger(this.getClass().getName()).fine("This molecule is considered NOT complex");
					CIPToolMod.label(cimp.getContainer());
				}

               for (int i = 0; i < container.getAtomCount(); i++) {
                   IAtom ai =container.getAtom(i);
                   IAtom ain =cimp.getContainer().getAtom(i);
                   Object p = ain.getProperty(CDKConstants.CIP_DESCRIPTOR);
                   
                   ai.removeProperty(CDKConstants.CIP_DESCRIPTOR);
                   if(p!=null) {
                	   ai.setProperty(CDKConstants.CIP_DESCRIPTOR, p);   
                   }
                   
               }
               
               for (int i = 0; i < container.getBondCount(); i++) {
                   IBond bi =container.getBond(i);
                   IBond bin =cimp.getContainer().getBond(i);
                   Object p = bin.getProperty(CDKConstants.CIP_DESCRIPTOR);
                   
                   bi.removeProperty(CDKConstants.CIP_DESCRIPTOR);
                   if(p!=null) {
                	   //For atropisomers
                	   bi.atoms().forEach(aa->aa.setProperty(CDKConstants.CIP_DESCRIPTOR, p));
                	   bi.setProperty(CDKConstants.CIP_DESCRIPTOR, p);   
                   }
                   
               }
               
               return 1;
               
            });
//            CIPTool.label(container);
//	    	container;
	    	
	    	this.doWithQueryFixes(()->{
	    	Stereocenters  centers   = Stereocenters.of(container);
	    	List<Integer> potentialSet = new ArrayList<>();
	    	List<Integer> undefinedSet = new ArrayList<>();
	    	for (int i = 0; i < container.getAtomCount(); i++) {
	    		switch(centers.stereocenterType(i)){
					case Non:
						break;
					case Para:
//						System.out.println("H");
//						break;
					case Potential:
						IAtom ai3=container.getAtom(i);
						// This is a weird kind of case where it really depends
						// on the circumstance. In general we care about any case where
						// if a wedge/hash were added this would change the
						// meaning of the molecule to something distinct.
						//
						// However this isn't perfect because sometimes adding a wedge or
						// dash in isolation won't change the meaning, but adding 2 or 3
						// will. For example 1,4-dimethylcyclohexane, or certain other
						// meso compounds. For those, we would be better off capturing
						// all potential as of yet undefined stereo centers and enumerating
						// the possibilities. If at least this specific center has more than one
						// possible R/S or r/s result from that enumerated set, then it's a
						// real stereocenter. If not, it's not a real stereocenter.
						
						if(ai3.getProperty(CDKConstants.CIP_DESCRIPTOR) == null){
							//The problem here is, if there IS a defined wedge/dash
							//it shouldn't be considered a potential center
							CdkAtom cat=new CdkAtom(container.getAtom(i),this);
							
							boolean hasStereo = cat.getBonds().stream()
					        .filter(b->b.getBondType().equals(BondType.SINGLE))
					        .filter(b->!b.getStereo().equals(Bond.Stereo.NONE))
					        .findAny().isPresent();
							if(!hasStereo) {
								potentialSet.add(i);
								undefinedSet.add(i);
								if(!deepChirality) {
									container.getAtom(i).setProperty(CDKConstants.CIP_DESCRIPTOR, "EITHER");
								}
							}
						}
						
						break;
					case True:
						IAtom ai=container.getAtom(i);
						if(ai.getProperty(CDKConstants.CIP_DESCRIPTOR) == null){
							// we only really set it if it's tetrahedral CARBON
							// right now,
							boolean bailout=false;
							//we have to go through the container because some atom implementations
                            //throw unsupport operation exceptions if we ask them for their bonds
							for(IBond ib:container.getConnectedBondsList(ai)){
								if(!Order.SINGLE.equals(ib.getOrder()) 
										
//								&& !Order.DOUBLE.equals(ib.getOrder())
										){
									bailout=true;
									break;
								}								
							}
							
							if(bailout){
								// axial cases are a little weird here
								// in general CDK marks all SP2-hybrid entries in rings
								// of < 9 as true stereocenters. I don't know why.
								
//								System.out.println("It's real but ...");
//								System.out.println(centers.elementType(i) + ":" + centers.stereocenterType(i) + ":" + ai.getSymbol());
								
								continue;
							}
							if("N".equals(ai.getSymbol())){
								if(ai.getFormalCharge()==0){
									//only charged nitrogens allowed
									continue;
								}
							}
							
							container.getAtom(i).setProperty(CDKConstants.CIP_DESCRIPTOR, "EITHER");
							undefinedSet.add(i);
						}
						break;
					default:
						break;
	    		}
	    	}

				Logger.getLogger(this.getClass().getName()).fine(String.format("number of potential centers: %d; number of undefined center: %d\n",
						potentialSet.size(), undefinedSet.size()));
	    	//TODO fix for potential cases
//	    	container.getAtom(i);
	    	
//	    	BitSet.valueOf(null);
	    	
	    	
	    	if(deepChirality && !potentialSet.isEmpty() && undefinedSet.size() <= maxUndefinedStereoCenters) {
	    		CdkChemicalImpl cimp2=this.deepCopy();
	    		cimp2.setDeepChirality(false);
				Chemical c22 = new Chemical(cimp2);
		    	
		    	Set<Integer> isDefinable = new HashSet<>();
		    	for(long ii=0;ii<Math.pow(2, undefinedSet.size());ii++) {
		    		cimp2.setDirty();
					BitSet bs= BitSet.valueOf(new long[] {ii});
		    		for(int ks=0;ks<undefinedSet.size();ks++) {
	//	    			int rs = Integer.
		    			if(bs.get(ks)) {
		    				c22.getAtom(undefinedSet.get(ks))
		    				.setChirality(Chirality.R);
		    			}else {
		    				c22.getAtom(undefinedSet.get(ks))
		    				.setChirality(Chirality.S);
		    			}	    			
		    		}
		    		for(int ks=0;ks<undefinedSet.size();ks++) {
		    			Chirality c=c22.getAtom(undefinedSet.get(ks)).getChirality();
//		    			System.out.print(c.toString());
		    			if(c.isDefined()) {
		    				isDefinable.add(undefinedSet.get(ks));
		    			}
		    		}
//		    		System.out.println();
		    		if(isDefinable.size()==undefinedSet.size()) {
		    			break;
		    		}
		    	}
		    	for(int pi:potentialSet) {
		    		if(isDefinable.contains(pi)) {
		    			container.getAtom(pi).setProperty(CDKConstants.CIP_DESCRIPTOR, "EITHER");	
		    		}
		    	}	    
	    	} else {
				for(int pi:potentialSet) {
					container.getAtom(pi).setProperty(CDKConstants.CIP_DESCRIPTOR, "EITHER");
				}
			}
	    	
	    	},false);
	    	
	    	return null;
    }
    	);
    CachedSupplier<Integer> ringsSearcherSupplier = CachedSupplier.of(()->{
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
					return 1;
				} catch (CDKException | IllegalArgumentException e) {
					Logger.getLogger(this.getClass().getName()).fine(String.format("Error processing rings in molecule with formula %s (%s; ring size limit %d)",
							this.getFormula(), e.getMessage(), Math.min(container.getAtomCount(),12)));
					return 0;
				}
		    	

		
    });


	private final CachedSupplierGroup cachedSupplierGroup = new CachedSupplierGroup();
	private final ChemicalSource source;
	
	
	public CdkChemicalImpl(IAtomContainer container, Supplier<? extends ChemicalSource> source) {
		this(container, source.get());
	}
	public CdkChemicalImpl(IAtomContainer container, ChemicalSource source) {

//		for(int i=0; i< container.getAtomCount(); i++){
//			IAtom a = container.getAtom(i);
//			if(a.getSymbol().startsWith("_R")){
//				container.setAtom(i, new PseudoAtom(a));
//			}
//		}
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
		cachedSupplierGroup.add(complexitySupplier);

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
		    if(a==null)continue;
		    a=getAtom(a.getAtomIndexInParent());
//		    System.out.println("A:"+a);
			for(Bond b : a.getBonds()) {
				gov.nih.ncats.molwitch.Bond.Stereo oldStereo = b.getStereo();
				gov.nih.ncats.molwitch.Bond.Stereo newStereo = oldStereo.flip();
				
				if(oldStereo !=newStereo) {
//					IAtom center = CdkAtom.getIAtomFor(s.getCenterAtom());
					b.setStereo(newStereo);
				}
			}
		}
	}

	public void setDeepChirality(boolean chir) {
		this.deepChirality=chir;
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
	    // A note here: this method used to rely on a deep "clone" of
	    // CDK's AtomContainer. However it turns out that clone isn't
	    // as deep as CDK suggests it is, and mutations done to the cloned
	    // container are often also done to the original container.
	    //
	    // Instead, this makes a deep clone by copying to molfile and importing
	    // again. This is not ideal, obviously. Furthermore, we can't just
	    // export to SD format since information regarding aromaticity is lost
	    // on toSD sometimes (it's not lost on toMol).
	    //
	    // TODO: debug CDK directly to fix
	    //
	    
//		try {
//		    Chemical wc = (new Chemical(this));
//		    String b = wc.toMol();
//		    Chemical p =Chemical.parse(b);
//		    wc.getPropertyIterator().forEachRemaining(ew->{
//		        if(ew.getValue()!=null) {
//		            p.setProperty(ew.getKey(), ew.getValue());
//		        }
//		    });
//
//            CdkChemicalImpl imp=(CdkChemicalImpl) p.getImpl();
//            //It's unclear if this is necessary
//		    for(int i=0;i<this.container.getAtomCount();i++) {
//		        imp.container.getAtom(i).setImplicitHydrogenCount(this.container.getAtom(i).getImplicitHydrogenCount());
//		    }
//			return new CdkChemicalImpl(imp.container,source);
//		} catch (Exception e) {
//			throw new IllegalStateException(e);
//		}
	    
        try {
            return new CdkChemicalImpl(container.clone(), source);
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException(e);
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

	protected void setDirty() {

		cachedSupplierGroup.resetCache();
		container.notifyChanged();
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

	    // This is a little slap-dash. Forcing hydrogens to be explicit
	    // requires that they had previously been set to be implicit somehow.
	    // MolWitch CDK has a complex and inconsistent relationship with 
	    // implicit hydrogen counts. This case below forces some preliminary
	    // detection flags, then some forced implicit H count settings,
	    // and THEN attempts to turn them into explicit H.
	    //
	    // TODO: It is likely that not all these steps are necessary every
	    // time, and more research is needed to be consistent.

        setImplicitHydrogens();

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
	    // As with the code that makes H's EXPLICIT, this code
        // also starts by calculating implicit counts. When this is necessary
	    // isn't clear.
		try {
		    setImplicitHydrogens();
		    try {
		    	AtomContainerManipulator.suppressHydrogens(container);
		    }catch(Exception e) {
		    	// suppress for now, it tends to fail if there is a query atom
		    	// but implicit Hs mean less there anyway
		    }
		}catch(Exception e) {
			
		    throw new RuntimeException(e);
		}
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
		for(IBond bond : container.getConnectedBondsList(atom)){
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
		double d=Double.NaN;
		try {
			d = AtomContainerManipulator.getMass(container, AtomContainerManipulator.MolWeight);
		}catch(Exception e){

		}
		if(Double.isNaN(d) || (d<0.01 && this.getAtomCount()>0)) {

			double off = 0;

			Map<IAtom, Integer> realAtomicNums = new HashMap<IAtom, Integer>();


			for (IAtom a : container.atoms()) {
				Integer an = a.getAtomicNumber();
				if (an == null) {
					a.setAtomicNumber(2);
					realAtomicNums.put(a, an);
					continue;
				}
				if (cdkMissing.get(an)) {
					double m = NISTIsotopeFactory.INSTANCE.getMostAbundant(an)
							.filter(mm -> mm.getIsotopicComposition() != null)
							.map(aa -> aa.getRelativeAtomicMass().getValue().doubleValue())
							.orElseGet(() -> {

								if (mostStable[an] != 0) return (double) mostStable[an];
								List<Isotope> ilist = NISTIsotopeFactory.INSTANCE.getIsotopesFor(an).stream()
										.sorted(Comparator.comparing(an1 -> an1.getRelativeAtomicMass().getValue().doubleValue()))
										.collect(Collectors.toList());
								;
								if (ilist.isEmpty()) return 0.0;
								Isotope iso = ilist.get(ilist.size() / 2);

								return (double) Math.round(iso.getRelativeAtomicMass().getValue().doubleValue());
							});

					off += m;
					a.setSymbol("R");
					realAtomicNums.put(a, an);
				}
			}
			try {
				double base = AtomContainerManipulator.getMass(container, AtomContainerManipulator.MolWeight);


				d = base + off;

				//this means there's an isotope issue
				//the fix, for now, is this ugly one

			} finally {
				realAtomicNums.forEach((a, i) -> {
					a.setAtomicNumber(i);
				});
			}
		}
		return d;
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
	        	 //sometimes ruins other calculations
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
	public void generateCoordinates() throws MolwitchException{
		try {
			StructureDiagramGenerator coordinateGenerator = new StructureDiagramGenerator(container);
			doWithQueryFixes(coordinateGenerator::generateCoordinates,false);
			container = coordinateGenerator.getMolecule();
		}catch(Exception e) {
			throw new MolwitchException(e.getMessage(), e);
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
			    makeStereoElms();

			}

		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	private void makeStereoElms() {
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
				
			    // This turns out to have some significant issues, and doesn't
			    // properly kekulize many 5-membered rings when
			    // implicit/explicit H counts are lost (as is often the case in
			    // transit). Need a simple way to recover those.
			    
				Kekulization.kekulize(container);

				//kekulize doesn't touch the aromatic bond flags
				//so we want to I guess?
				for(IBond bond : container.bonds()){
					bond.setIsAromatic(false);
				}
				isAromatic=false;
			},false);
			
		} catch (Exception e) {
			// TODO: probably need a DEBUG flag to suppress this / decide to throw
			// It gets thrown/printed a lot.
			
//			e.printStackTrace();
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
			if(CdkUtil.isPseudoAtom(iAtom)){
				continue;
			}
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

    public <T> T withModifiedForm(Function<Chemical, T> calc) {
    	
    	
        Chemical c = new Chemical(this);

//        if(true)return calc.apply(c);
        
        List<Tuple<Runnable,Runnable>> operations = new ArrayList<>();


        c.atoms()
        .filter(at->!at.isPseudoAtom() && !at.isQueryAtom() && ! at.isRGroupAtom())
        .filter(at->at.getAtomicNumber()>=11)
        .map(at->Tuple.of(at,at.getBonds().stream()
                .filter(bb->bb.getBondType().getOrder()==2)
                .filter(bb->{
                    Atom oat=bb.getOtherAtom(at);
                    if(oat.getSymbol().equals("O")) {
                        if(oat.getBondCount()==1)return true;
                    }
                    return false;
                } )
                .collect(Collectors.toList())))
        .filter(t->t.v().size()>0)
        .forEach(t->{
            Atom at = t.k();

            t.v().forEach(b->{
                BondType obt=b.getBondType();
                Atom oat = b.getOtherAtom(at);

                operations.add(Tuple.of(()->{
                    //doer
                    b.setBondType(BondType.SINGLE);
                    at.setCharge(at.getCharge()+1);
                    oat.setCharge(oat.getCharge()-1);
                },()->{
                    //undoer
                    b.setBondType(obt);
                    at.setCharge(at.getCharge()-1);
                    oat.setCharge(oat.getCharge()+1);
                }));
            });

            at.getBonds().stream().forEach(b->{
                if(b.getBondType().getOrder()==1) {
                    Atom oat=b.getOtherAtom(at);
                    if(oat.getSymbol().equals("O")) {
                        int bcount = oat.getBondCount();

                        //implicit H
                        if(bcount==1) {
                            if(oat.getCharge()==0) {
                                operations.add(Tuple.of(()->{
                                    //doer
                                    oat.setCharge(oat.getCharge()-1);

                                    at.setCharge(at.getCharge()+1);
                                },()->{
                                    //undoer                              
                                    at.setCharge(at.getCharge()-1);
                                    oat.setCharge(oat.getCharge()+1);
                                }));
                            }else {
                                if(oat.getCharge()==-1) {
                                    operations.add(Tuple.of(()->{
                                        //doer
                                        //                                      oat.setCharge(oat.getCharge()-1);

                                        at.setCharge(at.getCharge()+1);
                                    },()->{
                                        //undoer                              
                                        at.setCharge(at.getCharge()-1);
                                        //                                      oat.setCharge(oat.getCharge()+1);
                                    }));
                                }
                            }
                        }else if ( bcount==2 ) {

                            Atom hat = oat.getNeighbors().stream().filter(aa->aa.getSymbol().equals("H")).findAny().orElse(null);

                            //explicit H
                            if(hat!=null) {
                                Optional<Bond> existBond =(Optional<Bond>) hat.bondTo(oat);
                                operations.add(Tuple.of(()->{
                                    //doer
                                    Bond ob = existBond.get();
                                    c.removeBond(ob);
//                                    c.removeAtom(hat);
                                    oat.setCharge(oat.getCharge()-1);
                                    at.setCharge(at.getCharge()+1);
                                },()->{
                                    //undoer                              
                                    at.setCharge(at.getCharge()-1);
                                    oat.setCharge(oat.getCharge()+1);
                                    Bond ob = existBond.get();
//                                    c.addAtom(hat);
                                    c.addBond(ob);

                                }));
                            }
                        }
                    }
                }
            });

        });

        try {
            for(int i=0;i<operations.size();i++) {
                operations.get(i).k().run();
            }       

            try {
//                Chemical useC = Chemical.parse(c.toSd());
                  Chemical useC = c;
                return calc.apply(useC);
            } catch (Exception e) {
                return calc.apply(c);
                // TODO Auto-generated catch block
            }
        }finally {
            //Undo all of the modifications
            for(int i=operations.size()-1;i>=0;i--) {
                operations.get(i).v().run();
            }
        }

    }

    @Override
    public List<TetrahedralChirality> getTetrahedrals() {
        return withModifiedForm(cc->((CdkChemicalImpl)cc.getImpl()).getTetrahedrals1());
    }
    
	public List<TetrahedralChirality> getTetrahedrals1() {

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
			case r:
			case s:
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
		setImplicitHydrogens(((CdkAtom)b.getAtom1()).getAtom());
        setImplicitHydrogens(((CdkAtom)b.getAtom2()).getAtom());
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
			case AROMATIC : order = Order.UNSET; break; 
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
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
			hydrogenAdder.addImplicitHydrogens(container);
		} catch (CDKException e) {
//			e.printStackTrace();
		    throw new RuntimeException(e);
		}
	}
	protected int setImplicitHydrogens(IAtom atom){
	    return setImplicitHydrogens(getCdkAtomFor(atom));
	}
	protected int setImplicitHydrogens(CdkAtom catom){

        IAtom atom = catom.getAtom();
	    try {
	        atom.setImplicitHydrogenCount(null);
	        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
            hydrogenAdder.addImplicitHydrogens(container, atom);
            Integer cc= atom.getImplicitHydrogenCount();
            
            int tb = catom.getBondCount();
            
            int aromSNG=0;
            int aromDBL=0;
            int aromUNSET=0;
            
            for(Bond bbb:catom.getBonds()) {
                CdkBond cb = (CdkBond)bbb;
                IBond bb=cb.getBond();
                if(bb.isAromatic()) {
                    if(bb.getOrder()==Order.UNSET) {
                        aromUNSET++;
                    }else if(bb.getOrder()==Order.SINGLE) {
                        aromSNG++;
                    }else if(bb.getOrder()==Order.DOUBLE) {
                        aromDBL++;
                    }
                }
            }
            if((aromUNSET==3 && tb==3) || 
               (aromUNSET==2 && tb==3)){
                atom.setImplicitHydrogenCount(0);
                return 0;
            }
            if(aromUNSET==2 && tb==2) {
                if(catom.getSymbol().equals("C")) {
                    atom.setImplicitHydrogenCount(1);
                    return 1;
                }else {
                    atom.setImplicitHydrogenCount(0);
                }
            }
            
            return cc;
        } catch (CDKException e) {
            
//            return atom.getImplicitHydrogenCount();
            throw new RuntimeException("error computing implicit H count for atom " +atom, e);
        }
	    
	    //compute it
	    
	    
//		try {
//		    
//		    //TODO: This pipeline needs to be rewritten. It doesn't make sense.
//		    // ...? This logic doesn't make sense.
//		    
//			CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(CHEM_OBJECT_BUILDER);
//			 IAtomType type = matcher.findMatchingAtomType(container, atom);
//			
//			 AtomTypeManipulator.configure(atom, type);
//			//includes explicit and implicit Hs
//			 Integer neighborCount = atom.getFormalNeighbourCount();
//			 if(neighborCount ==null){
//				 neighborCount = 0;
//			 }
//			 double implicitH = neighborCount;
////			 System.out.println(implicitH + "  " + atom);
//			 for(IBond bond : container.getConnectedBondsList(atom)){
////			 	System.out.println("\t"+bond.getOrder());
//			     if(bond.isAromatic()) {
//			         implicitH-=1.5;
//			     }else {
//    				 switch(bond.getOrder()){
//    				 	case SINGLE : implicitH--; break;
//    				 	case DOUBLE : implicitH-=2; break;
//    				 	case TRIPLE : implicitH-=3; break;
//    				 	
//    				 	//???
//    				 	//WHY?
//    				 	case QUADRUPLE :implicitH-=4; break;  
//    				 	case QUINTUPLE :implicitH-=5; break;
//    				 	case SEXTUPLE :implicitH-=6; break;
//    				 	default : // ? TODO do we subtract 1?
//    				 }
//			     }
//			 }
////			 System.out.println(implicitH);
//			atom.setImplicitHydrogenCount(Math.max(0, (int)implicitH));
//			return Math.max(0, (int)implicitH);
//			 
//		}catch(CDKException e){
//			throw new RuntimeException("error computing implicit H count for atom " +atom, e);
//		}
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
			
			return Chirality.R.equals(c) || Chirality.S.equals(c) ||
					Chirality.r.equals(c) || Chirality.s.equals(c)
					;
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

	@Override
	public List<String> getSGroupWarnings() {
		List<String> messages = new ArrayList<>();
		AtomicInteger counter = new AtomicInteger(0);
		this.getSGroups().forEach(sg->{
			Assert.assertEquals(2, sg.getBrackets().size());
			Logger.getLogger(CdkChemicalImpl.class.getName()).info(String.format("looking at SGroup %d with %d atoms bracket 0 x: %.2f; y: %.2f; bracket 1 x: %.2f; y: %.2f\n",
					counter.incrementAndGet(), sg.getAtoms().count(), sg.getBrackets().get(0).getPoint1().getX(), sg.getBrackets().get(0).getPoint1().getY(),
					 sg.getBrackets().get(1).getPoint1().getX(), sg.getBrackets().get(1).getPoint1().getY()));
			double lowerX = sg.getBrackets().get(0).getPoint1().getX();
			double upperX = sg.getBrackets().get(1).getPoint1().getX();
			double lowerY = sg.getBrackets().get(0).getPoint1().getY();
			double upperY = sg.getBrackets().get(1).getPoint1().getY();
			for(int i = 0; i < this.getAtomCount(); i++) {
				Logger.getLogger(CdkChemicalImpl.class.getName()).info(String.format("	atom %d", i));
				Atom atom = this.getAtom(i);
				if( !sg.getAtoms().anyMatch( a-> a.equals(atom))){
					Logger.getLogger(CdkChemicalImpl.class.getName()).info(String.format("atom %s %d\n",
							atom.getSymbol(), i));
					if((atom.getAtomCoordinates().getX() >= lowerX && atom.getAtomCoordinates().getX() <= upperX)
						&& (atom.getAtomCoordinates().getY() >= lowerY && atom.getAtomCoordinates().getY() <= upperY)){
						String message = String.format("atom %s (%d) is within the bounds of an SGroup but not part of the SGroup",
								atom.getSymbol(), i);
						Logger.getLogger(CdkChemicalImpl.class.getName()).warning(message);
						messages.add(message);
					}
				}
			}
		});
		return messages;
	}
}
