/***************************************************************************
  * Title:          
  * Author:         Kyowon Jeong
  * Last modified:  
  *
  * Copyright (c) 2008-2009 The Regents of the University of California
  * All Rights Reserved
  * See file LICENSE for details.
  ***************************************************************************/
package uninovo;

import java.util.ArrayList;
import java.util.HashMap;

import uninovo.util.Composition;
import uninovo.util.IonType;
import uninovo.util.Peak;
import uninovo.util.Peptide;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;

/**
 * The Class PeakGenerator generates theoretical peaks of a (fragmented) peptide and an ion.
 */
public class PeakGenerator {
	
	/** The Constant VirtualPeakIntensity. */
	static final float VirtualPeakIntensity = Float.POSITIVE_INFINITY;
	
	/**
	 * Gets the charge changed peak.
	 *
	 * @param p the peak
	 * @param chargeBefore the charge before
	 * @param chargeOffset the charge offset
	 * @return the charge changed peak
	 */
	static public Peak getChargeChangedPeak(Peak p, int chargeBefore, int chargeOffset){
		float newmz;
		int chargeAfter = chargeOffset + chargeBefore;
		
		if(chargeBefore == chargeAfter){
			newmz = p.getMz();
		}else{
			float mass = (float) ((p.getMz() - Composition.PROTON) * chargeBefore);
			newmz = (float) (mass/chargeAfter + Composition.PROTON);
		}
		
		return new Peak(newmz, VirtualPeakIntensity, chargeAfter);
	}
	
	/**
	 * Gets the complementary peak.
	 *
	 * @param p the peak
	 * @param charge the charge
	 * @param pm the pm
	 * @return the complementary peak
	 */
	static public Peak getComplementaryPeak(Peak p, int charge, float pm){
		return new Peak((float) (pm/charge - p.getMz() + 2 * Composition.PROTON), VirtualPeakIntensity, charge);
	}
	
	/**
	 * Gets the complementary peak.
	 *
	 * @param p the peak
	 * @param charge the charge
	 * @param spec the spec
	 * @return the complementary peak
	 */
	static public Peak getComplementaryPeak(Peak p, int charge, Spectrum spec){
		return new Peak((float) (spec.getParentMass()/charge - p.getMz() + 2 * Composition.PROTON), VirtualPeakIntensity, charge);
	}
	
	/**
	 * Gets the prefix mass.
	 *
	 * @param p the peak
	 * @param ion the ion
	 * @param spec the spec
	 * @return the prefix mass
	 */
	static public float getPrefixMass(Peak p, IonType ion, Spectrum spec){
		float mass = ion.getMass(p.getMz());
		
		if(ion instanceof IonType.SuffixIon){
			mass = spec.getPeptideMass() - mass;
		}
		
		return mass;
	}
	
	/** The ion theoretical base peak map. */
	private HashMap<IonType, ArrayList<Peak>> ionTheoreticalBasePeakMap = null; // theoretical base peaks
	
	/** The peptide. */
	private Peptide pep;
	
	/**
	 * Instantiates a new peak generator.
	 *
	 * @param spec the spectrum
	 */
	public PeakGenerator(Spectrum spec){
		pep = spec.getAnnotation();
	}
	
	/**
	 * Gets the theoretical complementary peak.
	 *
	 * @param p the peak
	 * @param charge the charge
	 * @return the theoretical complementary peak
	 */
	Peak getTheoreticalComplementaryPeak(Peak p, int charge){
		if(pep == null) return null;
		return new Peak((float) (pep.getParentMass()/charge - p.getMz() + 2 * Composition.PROTON), VirtualPeakIntensity, charge);
	}

	/**
	 * Gets the theoretical peaks.
	 *
	 * @param ion the ion
	 * @return the theoretical peaks
	 */
	ArrayList<Peak> getTheoreticalPeaks(IonType ion){
		ArrayList<Peak> peaks;
		
		if(pep == null) return null;
		
		if(ionTheoreticalBasePeakMap == null)
			ionTheoreticalBasePeakMap = new HashMap<IonType, ArrayList<Peak>>();
		
		if((peaks = ionTheoreticalBasePeakMap.get(ion)) != null)
			return peaks;
		else peaks = new ArrayList<Peak>();
		
		if(ion instanceof IonType.PrefixIon){
			peaks = this.getTheoreticalPrefixPeaks(ion.getCharge(), ion.getOffset() * ion.getCharge());
		}else if(ion instanceof IonType.SuffixIon){
			peaks = this.getTheoreticalSuffixBasePeaks(ion.getCharge(), ion.getOffset() * ion.getCharge());
		}else if(ion instanceof IonType.PrecursorIon){
			peaks.add(this.getTheoreticalPrecursorPeak(ion.getCharge(), ion.getOffset() * ion.getCharge()));
		}
		
		ionTheoreticalBasePeakMap.put(ion, peaks);
		
		return peaks;
	}
	
	/**
	 * Gets the theoretical precursor peak.
	 *
	 * @param charge the charge
	 * @return the theoretical precursor peak
	 */
	Peak getTheoreticalPrecursorPeak(int charge){ return getTheoreticalPrecursorPeak(charge, 0); }

	/**
	 * Gets the theoretical precursor peak.
	 *
	 * @param charge the charge
	 * @param massOffset the mass offset
	 * @return the theoretical precursor peak
	 */
	Peak getTheoreticalPrecursorPeak(int charge, float massOffset) { 
		if(pep == null) return null;
		return new Peak((float) ((pep.getParentMass() + massOffset)/charge + Composition.PROTON), VirtualPeakIntensity, charge); 
	}
	
	/**
	 * Gets the theoretical prefix peaks.
	 *
	 * @param charge the charge
	 * @return the theoretical prefix peaks
	 */
	ArrayList<Peak> getTheoreticalPrefixPeaks(int charge){ return getTheoreticalPrefixPeaks(charge, 0); }
	
	/**
	 * Gets the theoretical prefix peaks.
	 *
	 * @param charge the charge
	 * @param massOffset the mass offset
	 * @return the theoretical prefix peaks
	 */
	ArrayList<Peak> getTheoreticalPrefixPeaks(int charge, float massOffset){
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		
		if(pep == null) return null;
		
		for(float m : pep.getPRMMasses(true, 0)){
			peaks.add(new Peak((m + massOffset)/charge, VirtualPeakIntensity, charge));
		}
		
		return peaks;
		
	}
	
	/**
	 * Gets the theoretical suffix peaks.
	 *
	 * @param charge the charge
	 * @param massOffset the mass offset
	 * @return the theoretical suffix peaks
	 */
	ArrayList<Peak> getTheoreticalSuffixBasePeaks(int charge, float massOffset){
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		
		if(pep == null) return null;
		
		for(float m : pep.getPRMMasses(false, 0)){
			peaks.add(new Peak((m + massOffset)/charge, VirtualPeakIntensity, charge));
		}
		
		return peaks;
		
	}
	
	/**
	 * Gets the theoretical suffix peaks.
	 *
	 * @param charge the charge
	 * @return the theoretical suffix peaks
	 */
	ArrayList<Peak> getTheoreticalSuffixPeaks(int charge){ return getTheoreticalSuffixBasePeaks(charge, 0); }
	
	/**
	 * Checks if a peak is explained by an ion.
	 *
	 * @param p the peak
	 * @param ion the ion
	 * @param tol the MS2 tol
	 * @param pmtol the MS1 tol
	 * @return true, if is explained by
	 */
	public boolean isExplainedBy(Peak p, IonType ion, Tolerance tol, Tolerance pmtol){
		boolean isPrecursorIon = (ion instanceof IonType.PrecursorIon);
		Tolerance tmptol = isPrecursorIon? pmtol : tol;
		for(Peak peak : getTheoreticalPeaks(ion)){
			float diff = Math.abs(p.getMz() - peak.getMz());
			if(diff <= tmptol.getToleranceAsDa(p.getMz() * ion.getCharge())){
				return true;
			}
		}
		return false;
	}
}
