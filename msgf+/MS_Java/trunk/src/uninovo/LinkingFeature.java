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

import uninovo.util.AminoAcid;
import uninovo.util.AminoAcidSet;
import uninovo.util.Peak;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
 
/**
 * The Class LinkingFeature.
 */
public class LinkingFeature extends Feature{
	
	/** The current iteration number. */
	static private int currentIterationNum = -1;
	
	/** The current spectrum - just for speed-up. */
	static private Spectrum currentSpec = null;
	
	/**
	 * Parses the file string.
	 *
	 * @param s the s
	 * @param aaSet the aa set
	 * @return the linking feature
	 */
	static public LinkingFeature parseFileString(String s, AminoAcidSet aaSet){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseSpectrumParameter(token[2]);
		PeakParameter ppar = PeakParameter.parsePeakParameter(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		int peakIntensityRatio = Integer.parseInt(token[5]);
		boolean isPresent = token[6].equals("p");
		
		return new LinkingFeature(spar, ppar, basePeakCharge, isPresent, iterationNum, peakIntensityRatio, aaSet);
	}
	
	/** The amino acid set. */
	private AminoAcidSet aaSet;

	/** The current hold. */
	private boolean currentHold = true;
	
	/** The current peak. */
	private Peak currentPeak = null;
	/**
	 * Instantiates a new linking feature.
	 *
	 * @param spar the spectrum parameter
	 * @param ppar the peak parameter
	 * @param basePeakCharge the base peak charge
	 * @param isPresent the peak that make this feature satisfy is present
	 * @param iterationNum the iteration number
	 * @param peakIntensityRatio the peak intensity ratio
	 * @param aaSet the amino acid set
	 */
	public LinkingFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge, boolean isPresent,
			int iterationNum, int peakIntensityRatio, AminoAcidSet aaSet) {
		super(spar, ppar, basePeakCharge, peakIntensityRatio, isPresent, iterationNum);
		this.setPeakIntensityRatio(peakIntensityRatio);
		this.aaSet = aaSet;
	}

	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof LinkingFeature) ) return false;
    	LinkingFeature con = (LinkingFeature)o;
    	
    	return this.getPeakIntensityRatio() == con.getPeakIntensityRatio() && super.equals(con);
	}


	@Override
	public int hashCode() {
		return super.hashCode()	* (this.getPeakIntensityRatio() + 12773);
	}

	@Override
	public boolean isSatisfiedBy(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol, Tolerance pmtol, int iterationNum) {
		if(currentSpec!=null && spec.getParentMass() == currentSpec.getParentMass() && spec.equals(currentSpec) && bp == currentPeak){
			return currentHold;
		}
		
		boolean hold = false;

		if(super.isSatisfiedBy(bp, iterationNum, spar, ppar)){	
			if(!spec.equals(currentSpec) || currentIterationNum != iterationNum){
				currentSpec = spec; currentIterationNum = iterationNum;
			}

			boolean match = false;
			for(AminoAcid aa : aaSet){
				for(int i=0; i<2; i++){
					float mz = bp.getMz() + aa.getMass()/bp.getCharge() * (i == 0 ? 1 : -1);
					float t = tol.getToleranceAsDa(bp.getMz() * this.getBasePeakCharge());
					for(Peak cp : spec.getPeakListByMassRange(mz - t, mz + t)){
						int ratio = PeakParameter.getPeakIntensityRatioNum(bp, cp, spec);
						if(ratio == this.getPeakIntensityRatio()){
							match = true;							
							break;
						}
					}
					if(match) break;
				}
				if(match) break;
			}
		
			
			hold = match == this.isPresent();	
		}
		currentPeak = bp;
		currentHold = hold;
		return hold;
	}

	@Override
	public String toFileString() {
		return "B\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString()+"\t" +this.getBasePeakCharge() + "\t" + this.getPeakIntensityRatio() + "\t" + (this.isPresent()? "p" : "a");
	}
	
	@Override
	public String toString() {
		return "Bridging - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge() + " Ratio: " + this.getPeakIntensityRatio() + " Present: " + this.isPresent();
	}




}
