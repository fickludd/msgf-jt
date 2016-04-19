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

import uninovo.util.Peak;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;

/**
 * The Class OffsetFeature describes an offset feature
 */
public class OffsetFeature extends Feature{
	
	/**
	 * Parses the file string.
	 *
	 * @param s the s
	 * @return the offset feature
	 */
	static public OffsetFeature parseFileString(String s){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseSpectrumParameter(token[2]);
		PeakParameter ppar = PeakParameter.parsePeakParameter(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		int peakIntensityRatio = Integer.parseInt(token[5]);
		boolean isPresent = token[6].equals("p");
		RelationBetweenPeaks gof = RelationBetweenPeaks.parseFileString(token[7]);
		
		return new OffsetFeature(spar, ppar, basePeakCharge, isPresent, iterationNum, peakIntensityRatio, gof);
	}
	
	/** The current hold. */
	private boolean currentIsSatisfied = true;


	/** The current peak. */
	private Peak currentPeak = null;
	
	/** The current spec. */
	private Spectrum currentSpec = null;
	
	/** The gof. */
	private RelationBetweenPeaks relation; // redundant infomation - should be optimized..
	

	
	/**
	 * Instantiates a new offset feature.
	 *
	 * @param spar the spec parameter
	 * @param ppar the peak parameter
	 * @param basePeakCharge the base peak charge
	 * @param isPresent is a peak that make this satisfies present
	 * @param iterationNum the iteration number
	 * @param peakIntensityRatio the peak intensity ratio
	 * @param relation the peak relation of this feature
	 */
	public OffsetFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge, boolean isPresent,
			int iterationNum, int peakIntensityRatio, RelationBetweenPeaks relation) {
		super(spar, ppar, basePeakCharge, peakIntensityRatio, isPresent, iterationNum);
		assert(basePeakCharge == relation.getBaseCharge());
		this.setRelation(relation);
	}

	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof OffsetFeature) ) return false;
    	OffsetFeature con = (OffsetFeature)o;
    	
    	return this.getPeakIntensityRatio() == con.getPeakIntensityRatio() && super.equals(con) && this.getRelation().equals(con.getRelation());
	}

	/**
	 * Gets the charge offset.
	 *
	 * @return the charge offset
	 */
	public int getChargeOffset() {return getRelation().getChargeOffset();}

	/**
	 * Gets the relation.
	 *
	 * @return the relation
	 */
	public RelationBetweenPeaks getRelation() {
		return relation;
	}
	
	@Override
	public int hashCode() {
		return super.hashCode() * (this.getPeakIntensityRatio() + 48701) * this.getRelation().hashCode();
	}
	
	/**
	 * Checks if this is complementary.
	 *
	 * @return true, if this is complementary
	 */
	public boolean isComplementary() {return getRelation().isComplementary();}
	
	@Override
	public boolean isSatisfiedBy(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol, Tolerance pmtol, int iterationNum) {
		if(currentSpec!=null && spec.getParentMass() == currentSpec.getParentMass() && spec.equals(currentSpec) && bp == currentPeak){
			return currentIsSatisfied;
		}
		
		boolean hold = false;

		if(super.isSatisfiedBy(bp, iterationNum, spar, ppar)){
			boolean match = false;
		
			for(Peak cp : getRelation().getSupportingPeaks(bp, spec, tol, pmtol)){
				if(getPeakIntensityRatio() == PeakParameter.getPeakIntensityRatioNum(bp, cp, spec)){
					match = true;
					break;
				}
			}
			hold = match == this.isPresent();	
		}
		
		currentSpec = spec;
		currentIsSatisfied = hold;
		currentPeak = bp;
		
		return hold;
	}
	
	/**
	 * Sets the relation.
	 *
	 * @param relation the new relation
	 */
	public void setRelation(RelationBetweenPeaks relation) {
		this.relation = relation;
	}

	@Override
	public String toFileString() {
		return "G\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString() +"\t" +this.getBasePeakCharge() + "\t" + this.getPeakIntensityRatio() + "\t" + (this.isPresent()? "p" : "a") + "\t" + this.getRelation().toFileString();
	}

	@Override
	public String toString() {
		return "GOF - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge() + " Ratio: " + this.getPeakIntensityRatio() + " Present: " + this.isPresent() + " GOF: " + this.getRelation();
	}

	
}
