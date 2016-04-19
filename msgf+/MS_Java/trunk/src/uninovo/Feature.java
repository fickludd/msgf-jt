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

import uninovo.util.AminoAcid;
import uninovo.util.Composition;
import uninovo.util.IonType;
import uninovo.util.Peak;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;

/**
 * An abstract class representing a feature
 * @author kyowon
 */
public abstract class Feature implements Comparable<Feature>{
	
	/**
	 * Returns the peak ion probability distribution
	 * @param bp peak
	 * @param features features that the peak satisfy
	 * @param ions ions
	 * @param spec spectrum
	 * @return the peak ion probability distribution
	 */
	static public HashMap<IonType, Float> getPeakIonProbabilityDistribution(Peak bp, ArrayList<Feature> features, ArrayList<IonType> ions, Spectrum spec){
		HashMap<IonType, Float> ret = new HashMap<IonType, Float>();
		 
		for(IonType ion : ions){ 
			ret.put(ion, getUnNormalizedPeakIonProbability(ion, features, bp , spec));
		}
		
		if(!ions.contains(IonType.NOISE)) ret.put(IonType.NOISE, getUnNormalizedPeakIonProbability(IonType.NOISE, features,  bp , spec));
		
		float sum = 0;
		for(IonType i : ret.keySet()){
			sum+=ret.get(i);
		}
	
		if(sum == 0){
			sum = 1;
			if(ret.containsKey(IonType.NOISE)) ret.put(IonType.NOISE, 1f);
		}
		for(IonType i : ret.keySet()){
			ret.put(i, ret.get(i)/sum);
		}
		return ret;
	}
	
	/**
	 * Returns the peak ion probability (unnormalized)
	 * @param ion ion
	 * @param features features that the peak satisfy
	 * @param bp peak
	 * @param spec spectrum
	 * @return the peak ion probability 
	 */
	static private float getUnNormalizedPeakIonProbability(IonType ion, ArrayList<Feature> features, Peak bp, Spectrum spec){
		float p = 1;
		IntensityFeature nf = (IntensityFeature) features.get(0).intensityFeature;
		
		p = nf.getProbability(ion);
		
		for(Feature feature : features){
			if(feature instanceof IntensityFeature) continue;	
			float mult = 1;
			if(bp.getMz() < minMzFor(ion) || bp.getMz() > maxMzFor(ion, spec)) mult = 0;
			else mult = feature.getProbability(ion);
			
			p*=mult;
		}

		return p;
	}
	
	/**
	 * Returns the maximum mz value for an ion and a spectrum 
	 * @param ion ion
	 * @param spec spectrum
	 * @return the maximum mz
	 */
	static public float maxMzFor(IonType ion, Spectrum spec){
		if(ion.equals(IonType.NOISE)) return Float.MAX_VALUE;
		return ion.getMz((float) (spec.getParentMass() - Composition.H2O - AminoAcid.getStandardAminoAcid('G').getMass()));
	}
	
	/**
	 * Returns the minimum mz value for an ion and a spectrum.
	 *
	 * @param ion ion
	 * @return the minimum mz
	 */
	static public float minMzFor(IonType ion){
		if(ion.equals(IonType.NOISE)) return -Float.MAX_VALUE;
		return ion.getMz(AminoAcid.getStandardAminoAcid('G').getMass());
	}
	
	// peak charge
	/** The base peak charge. */
	private int basePeakCharge;
	
	// peak ion probability distribution
	/** The distribution. */
	private HashMap<IonType, Float> distribution = null;
	
	// intensity feature of this feature
	/** The intensity feature. */
	private Feature intensityFeature = null;
	
	// if the other peak that makes the base peak satisfy this feature present or not
	/** The is present. */
	private boolean isPresent;
	
	// iteration number
	/** The iteration number. */
	private int iterationNumber;
	
	// number of peaks satisfying this feature - used only for training
	/** The number of satisfied peaks. */
	private HashMap<IonType, Float> numberOfSatisfiedPeaks = null;

	// peak intensity ratio value
	/** The peak intensity ratio. */
	private int peakIntensityRatio;
	
	// peak parameter
	/** The ppar. */
	private PeakParameter ppar;
	
	// spectrum parameter 
	/** The spar. */
	private SpectrumParameter spar;

	/**
	 * Constructor.
	 *
	 * @param spar spectrum parameter
	 * @param ppar peak parameter
	 * @param basePeakCharge base peak charge
	 * @param peakIntensityRatio peak intensity ratio
	 * @param isPresent if the other peak that makes the base peak satisfy this feature present or not
	 * @param iterationNumber the iteration number
	 */
	protected Feature(SpectrumParameter spar, PeakParameter ppar, int basePeakCharge, int peakIntensityRatio, boolean isPresent, int iterationNumber){
		this.spar = spar;
		this.ppar = ppar;
		this.basePeakCharge = basePeakCharge;
		this.peakIntensityRatio = peakIntensityRatio;
		this.isPresent = isPresent;
		this.iterationNumber = iterationNumber;
		if(this.distribution == null) this.distribution = new HashMap<IonType, Float>();
	}

	/**
	 * Used when training to count the number of ions satisfying this feature 
	 * All ions including Noise are counted.
	 * @param ion ion
	 */	
	public void addIonCount(IonType ion){ 
		if(numberOfSatisfiedPeaks == null) numberOfSatisfiedPeaks = new HashMap<IonType, Float>();
		Float num = numberOfSatisfiedPeaks.get(ion);
		if(num == null) num = 0f;
		num++;
		numberOfSatisfiedPeaks.put(ion, num);
	}
	
	/**
	 * Used when training to calculate the ion probabilities of this feature.
	 * Called once after adding all ion counts in the training dataset
	 * @return the total number of peaks satisfying this feature
	 */	
	public float calculateProbabilities(){
		if(numberOfSatisfiedPeaks == null) return 0;
		float sum = 0;
		for(IonType ion : numberOfSatisfiedPeaks.keySet()){
			sum += numberOfSatisfiedPeaks.get(ion);
		}
		
		for(IonType ion : numberOfSatisfiedPeaks.keySet()){
			float p = 1;
			if(this instanceof IntensityFeature){
				p = numberOfSatisfiedPeaks.get(ion)/sum;
			}else{
				p = numberOfSatisfiedPeaks.get(ion)/this.intensityFeature.numberOfSatisfiedPeaks.get(ion);
			}
			if(!(ion instanceof IonType.PrecursorIon)) distribution.put(ion, p);
		}
		return sum;
	}
	
	/**
	 * comparator - compares the divergences
	 * overrides compareTo
	 * @param o other feature
	 * @returns comparison result
	 */	
	@Override
	public int compareTo(Feature o) {
		return new Float(this.getDivergence()).compareTo(new Float(o.getDivergence()));
	}

	/**
	 * overrides equals
	 * @param o object
	 * @return if this and o equal or not
	 */	
	@Override
	public boolean equals(Object o){
		if(this == o) return true;
    	if(!(o instanceof Feature)) return false;
    	Feature con = (Feature)o;
    	return con.isPresent == this.isPresent &&
    		con.iterationNumber == this.iterationNumber &&
	    	con.basePeakCharge == this.basePeakCharge &&
	    	con.spar.equals(this.spar) &&
	    	con.ppar.equals(this.ppar);
	}

	/**
	 * Getter.
	 *
	 * @return the base peak charge
	 * @returns basePeakCharge
	 */
	public int getBasePeakCharge() {return basePeakCharge;}
	
	/**
	 * Getter.
	 *
	 * @return the base peak parameter
	 * @returns peak parameter
	 */
	public PeakParameter getBasePeakParameter() {return ppar;}
	
	/**
	 * Returns the divergence of this feature.
	 *
	 * @return the divergence
	 * @returns divergence
	 */
	public float getDivergence(){
		if(intensityFeature == null || this instanceof IntensityFeature) return 0;
		
		float kl = 0;
		
		for(IonType ion:this.distribution.keySet()){
			if(ion instanceof IonType.PrecursorIon) continue;
			
			float sum = 0;
			float nom = 0;
			for(IonType i : distribution.keySet()){
				float t = getProbability(i) * this.intensityFeature.getProbability(i);
				if(i.equals(ion)) nom = t;
				sum += t;
			}
			float p1 = nom/sum;
			
			
			//float p1 = getProbability(ion);
			if(p1 <= 0) continue;
			float p2 = intensityFeature.getProbability(ion);	
			
			assert(p2 !=0);
			
			kl+= (float)(p1 * Math.log(p1/p2) / Math.log(2) );
		}

		return kl;
	}
	
	/**
	 * Getter.
	 *
	 * @return the intensity feature
	 * @returns the intensity feature
	 */
	public Feature getIntensityFeature() {return intensityFeature;}
	
	/**
	 * Getter.
	 *
	 * @return the iteration num
	 * @returns iterationNum
	 */
	public int getIterationNum() {return iterationNumber;}
	
	/**
	 * Getter get peakIntensityRatio
	 * @return peakIntensityRatio 
	 */	
	public int getPeakIntensityRatio() {
		return peakIntensityRatio;
	}
	
	/**
	 * Returns the probability a peak satisfying this feature represents an ion (ion probability  of this feature)
	 * It is mu_f if this is not an intensityFeature and is gamma otherwise (in paper).
	 *
	 * @param ion ion
	 * @return the probability
	 * @returns the probability
	 */	
	public float getProbability(IonType ion){
		if(distribution.containsKey(ion))
			return distribution.get(ion);
		else return 0f;
	}
	
	/**
	 * Getter.
	 *
	 * @return the spectrum parameter
	 * @returns spectrum parameter
	 */
	public SpectrumParameter getSpectrumParameter() {return spar;}
	
	/**
	 * overrides hashCode
	 * @return hash code
	 */	
	@Override
	public int hashCode(){
		int t = this.iterationNumber + this.basePeakCharge << 3 + spar.hashCode() + ppar.hashCode();
		return this.isPresent? t : -t;
	}
	
	/**
	 * Getter.
	 *
	 * @return true, if is present
	 * @returns isPresent
	 */
	public boolean isPresent() {return isPresent;}
	
	/**
	 * Tests if this feature is satisfied by a peak
	 * @param bp a peak
	 * @param iterationNumber iterationNumber
	 * @param spar spectrum parameter
	 * @param ppar peak parameter
	 * @return if this is satisfied or not
	 */	
	protected boolean isSatisfiedBy(Peak bp, int iterationNumber, SpectrumParameter spar, PeakParameter ppar){	
		bp.setCharge(basePeakCharge);
		return iterationNumber == this.iterationNumber &&
		spar.equals(this.spar) &&
		ppar.equals(this.ppar);
	}
	
	/**
	 * An abstract method that tests if this feature is satisfied by a peak.
	 *
	 * @param bp a peak
	 * @param spec a spectrum
	 * @param spar spectrum parameter
	 * @param ppar peak parameter
	 * @param tol MS2 tolerance
	 * @param pmtol MS1 tolerance
	 * @param iterationNum the iteration num
	 * @return if this is satisfied or not
	 */	
	public abstract boolean isSatisfiedBy(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol, Tolerance pmtol, int iterationNum);
	
	/**
	 * Setter  
	 * @param feature an intensity feature
	 */
	public void setIntensityFeature(Feature feature){intensityFeature = feature;}
	
	/**
	 * Used when reading a parameter file to update the distribution of this feature
	 * distribution should include Noise probability
	 * @param distribution distribution
	 */	
	public void setIonProbMap(HashMap<IonType, Float> distribution){
		assert(distribution.containsKey(IonType.NOISE));
		this.distribution = distribution;
	}
	
	/**
	 * Setter set peakIntensityRatio.
	 *
	 * @param peakIntensityRatio the new peak intensity ratio
	 */	
	public void setPeakIntensityRatio(int peakIntensityRatio) {
		this.peakIntensityRatio = peakIntensityRatio;
	}
	
	/**
	 * an abstract method that returns a string representation of this for file output
	 * @return a string representation of this for file output
	 */	
	public abstract String toFileString();
	
	/**
	 * an abstract method overriding toString
	 * @return a string representation of this
	 */	
	@Override
	public abstract String toString();
	
}


