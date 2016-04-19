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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import uninovo.util.Peak;
import uninovo.util.Tolerance;

/**
 * The Class FeatureFrequencyFunction (FFF) shows the frequency of a feature
 */
public class FeatureFrequencyFunction {
	
	/**
	 * The Class FeatureFrequencyFunctionPeak.
	 */
	static class FeatureFrequencyFunctionPeak  implements Comparable<FeatureFrequencyFunctionPeak>{
    	
	    /** The fff. */
	    private RelationBetweenPeaks fff;
    	
	    /** The frequency. */
	    private float y;
    	
    	/**
	     * Instantiates a new FeatureFrequencyFunctionPeak.
	     *
	     * @param fff the fff
	     * @param y the frequency
	     */
	    FeatureFrequencyFunctionPeak(RelationBetweenPeaks fff, float y){
    		this.fff = fff; this.y = y;
    	}
    	
    	@Override
		public int compareTo(FeatureFrequencyFunctionPeak o) {
    		return new Float(this.y).compareTo(new Float(o.y));
    	}
    	
    	/**
	     * Gets the FFF
	     *
	     * @return the FFF
	     */
	    public RelationBetweenPeaks getFFF() { return fff; }
    	
	    /**
	     * Gets the frequency.
	     *
	     * @return the frequency
	     */
	    public float getFrequency() { return y; }
    }
	
	/** The Constant MAX - the maximum m/z offset. */
	static public final float MAX = 38;
	
	/** The Constant MIN. - the minimum m/z offset */
	static public final float MIN = -38;

	/**
	 * Gets the feature frequency function FFF - also can be seen as a OFF
	 *
	 * @param featureFrequencies the feature frequencies
	 * @param normalizer the normalizer
	 * @param threshold the threshold to select features
	 * @return the feature frequency function peaks
	 */
	static public ArrayList<FeatureFrequencyFunctionPeak> getFeatureFrequencyFunction(HashMap<RelationBetweenPeaks, Integer> featureFrequencies, float normalizer, float threshold){
		ArrayList<FeatureFrequencyFunctionPeak> fffPeaks = new  ArrayList<FeatureFrequencyFunctionPeak>();
		if(featureFrequencies == null || featureFrequencies.isEmpty()) return fffPeaks;
		
		if(threshold > 0 && threshold < 0.15f){
		
			ArrayList<Float> v = new ArrayList<Float>();
			for(RelationBetweenPeaks gof : featureFrequencies.keySet()){
				v.add((float)featureFrequencies.get(gof)/normalizer);			
				//sum+= (float)featureFrequencies.get(gof)/normalizer;
			}
			Collections.sort(v);
		
			threshold = Math.min(threshold, v.get(v.size()/2)*7);
		}
		
		ArrayList<FeatureFrequencyFunctionPeak> fffPeakstmp = new ArrayList<FeatureFrequencyFunctionPeak>();
		
		for(RelationBetweenPeaks feature : featureFrequencies.keySet()){
			float prob = (float)featureFrequencies.get(feature)/normalizer;
			
			if(prob >= threshold)
				fffPeakstmp.add(new FeatureFrequencyFunctionPeak(feature, prob));
		}
		
		Collections.sort(fffPeakstmp, Collections.reverseOrder());
		
		for(int i=0; i<Math.min(fffPeakstmp.size(), 100); i++){
			fffPeaks.add(fffPeakstmp.get(i));
		}				
		
		return fffPeaks;
	}
	
	
		
	/**
	 * Gets the features between bp and cps
	 *
	 * @param bp the bp
	 * @param cps the cps
	 * @param isComplementary the features are complementary
	 * @param chargeOffset the charge offset of the features
	 * @param tol the MS2 tol
	 * @return the features
	 */
	static public HashSet<RelationBetweenPeaks> getFeaturesBetween(Peak bp, ArrayList<Peak> cps, boolean isComplementary, int chargeOffset, Tolerance tol){
		HashSet<RelationBetweenPeaks> features = new HashSet<RelationBetweenPeaks>();
		for(Peak cp : cps){
			RelationBetweenPeaks feature = getRelationsBetween(bp, cp, isComplementary, chargeOffset, tol);
			if(feature != null) features.add(feature);
		}
		
		return features;
	}

	
	/**
	 * Gets the features between bp and cp
	 *
	 * @param bp the bp
	 * @param cp the cp
	 * @param isComplementary the features are complementary
	 * @param chargeOffset the charge offset of the features
	 * @param tol the MS2 tol
	 * @return the feature
	 */
	static public RelationBetweenPeaks getRelationsBetween(Peak bp, Peak cp, boolean isComplementary, int chargeOffset, Tolerance tol){
		// base peaks already have charge offset	
		float offset = (cp.getMz() - bp.getMz()) * bp.getCharge();
		if(offset > MAX || offset < MIN) return null;
		return new RelationBetweenPeaks(offset, isComplementary, chargeOffset, bp.getCharge() - chargeOffset, tol.getToleranceAsDa(500)*2);
	}
	
}
