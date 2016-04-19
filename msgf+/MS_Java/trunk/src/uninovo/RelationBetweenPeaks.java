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

import uninovo.util.Constants;
import uninovo.util.Peak;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;

/**
 * The Class RelationBetweenPeaks describes the relation between two peaks considering 
 * m/z offset charge offset complementary relation. This is used to define offset feature and for training.
 */
public class RelationBetweenPeaks {
	
	/**
	 * Parses the file string.
	 *
	 * @param s the s
	 * @return the RelationBetweenPeaks for training
	 */
	static public RelationBetweenPeaks parseFileString(String s){ 
		String[] token = s.split(" ");
		
		return new RelationBetweenPeaks(Float.parseFloat(token[0]), Integer.parseInt(token[1]) == 1,
				Integer.parseInt(token[2]), Integer.parseInt(token[3]));
	}
	
	/** The base charge. */
	private int baseCharge;
	
	/** The charge offset. */
	private int chargeOffset;
	
	/** is this complementary? */
	private boolean isComplementary;
	
	/** The offset. */
	private float offset;
	
	/**
	 * Instantiates a new relation between peaks.
	 *
	 * @param offset the offset
	 * @param isComplementary the is complementary
	 * @param chargeOffset the charge offset
	 * @param baseCharge the base charge
	 */
	private RelationBetweenPeaks(float offset, boolean isComplementary, int chargeOffset, int baseCharge){
		this.offset = offset;
		this.isComplementary = isComplementary;
		this.chargeOffset = chargeOffset;
		this.baseCharge = baseCharge;
	}
	
	/**
	 * Instantiates a new relation between peaks.
	 *
	 * @param offset the offset
	 * @param isComplementary the is complementary
	 * @param chargeOffset the charge offset
	 * @param baseCharge the base charge
	 * @param tolInDa the da tolerance of training data set
	 */
	public RelationBetweenPeaks(float offset, boolean isComplementary, int chargeOffset, int baseCharge, float tolInDa){
		if(tolInDa > 0){
			this.offset = (Math.round(offset*Constants.INTEGER_MASS_SCALER_HIGH_PRECISION))/Constants.INTEGER_MASS_SCALER_HIGH_PRECISION;		
		}
		else this.offset = offset;

		if(tolInDa >= 0.5) 
			this.offset = Math.round(offset*Constants.INTEGER_MASS_SCALER) / Constants.INTEGER_MASS_SCALER;
		
		this.isComplementary = isComplementary;
		this.chargeOffset = chargeOffset;
		this.baseCharge = baseCharge;
	}
	
	@Override
	public boolean equals(Object o){
    	if(this == o) return true;
    	if ( !(o instanceof RelationBetweenPeaks) ) return false;
    	RelationBetweenPeaks io = (RelationBetweenPeaks)o;
    	return io.offset == this.offset && io.baseCharge == this.baseCharge
    		&& io.isComplementary == this.isComplementary && io.chargeOffset == this.chargeOffset;
    }
	
	/**
	 * Gets the base charge.
	 *
	 * @return the base charge
	 */
	public int getBaseCharge() { return baseCharge; }
	
	/**
	 * Gets the charge offset.
	 *
	 * @return the charge offset
	 */
	public int getChargeOffset() { return chargeOffset; }
	
	/**
	 * Gets the offset.
	 *
	 * @return the offset
	 */
	public float getOffset(){ return offset; }
	
	/**
	 * Gets the peaks that make a base peak satisfy this peak.
	 *
	 * @param bp the base peak
	 * @param spec the spec
	 * @param tol the MS2 tol
	 * @param pmtol the MS1 tol
	 * @return the supporting peaks
	 */
	public ArrayList<Peak> getSupportingPeaks(Peak bp, Spectrum spec, Tolerance tol, Tolerance pmtol){
		if(bp.getCharge() != baseCharge) return new ArrayList<Peak>();
		
		Tolerance to = tol;
		
		if(isComplementary){
			if(tol.getToleranceAsDa(100) < pmtol.getToleranceAsDa(100)) to = pmtol;
		}
		
		float t = to.getToleranceAsDa(bp.getMz() * bp.getCharge());

		if(chargeOffset != 0)
			bp = PeakGenerator.getChargeChangedPeak(bp, bp.getCharge(), chargeOffset);
		if(isComplementary)
			bp = PeakGenerator.getComplementaryPeak(bp, bp.getCharge(), spec);
			
		float mz = bp.getMz() + offset / bp.getCharge();
		
		float s = to.getToleranceAsDa(mz * bp.getCharge());
		
		if(isComplementary && t < 0.25f){
			float  u = pmtol.getToleranceAsDa(mz * bp.getCharge());
			s = s>u? s : u;
			t += s;			
		}else	
			t = t>s? t : s;
		
		t = Math.min(t, 0.5f);
		ArrayList<Peak> ret = spec.getPeakListByMassRange(mz-t, mz + t);
		
		if(bp.getCharge() == spec.getCharge() && Math.abs(spec.getPrecursorPeak().getMz() - mz) < pmtol.getToleranceAsDa(mz))
			ret.add(spec.getPrecursorPeak());
	
		return ret;
	}
	
	@Override
	public int hashCode(){
		return this.baseCharge * (this.chargeOffset + 10) * (int)(this.offset) * (isComplementary? 1 : -1); 
	}
	
	/**
	 * Checks if is complementary.
	 *
	 * @return true, if is complementary
	 */
	public boolean isComplementary() { return isComplementary; }
	
	/**
	 * To file string.
	 *
	 * @return the string
	 */
	public String toFileString(){
		String s = this.offset + " " + (this.isComplementary? 1 : 0) + " " + this.chargeOffset + " " + this.baseCharge;
		return s;
	}
	
	@Override
	public String toString(){
		return "Off: " + this.offset + " Comp: " + this.isComplementary + " Charge: " + this.baseCharge + " ChargeOff: " + this.chargeOffset;
	}
	
}
