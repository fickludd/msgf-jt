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
 * The Class IntensityFeature describes an intensity feature
 */
public class IntensityFeature extends Feature{
	
	
	/**
	 * Parses the file string.
	 *
	 * @param s the s
	 * @return the intensity feature
	 */
	static public IntensityFeature parseFileString(String s){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseSpectrumParameter(token[2]);
		PeakParameter ppar = PeakParameter.parsePeakParameter(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		
		return new IntensityFeature(spar,	ppar, basePeakCharge, iterationNum);
	}

	/**
	 * Instantiates a new intensity feature.
	 *
	 * @param spar the spar
	 * @param ppar the ppar
	 * @param basePeakCharge the base peak charge
	 * @param iterationnum the iteration number
	 */
	public IntensityFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge,
			int iterationnum) {
		super(spar, ppar, basePeakCharge, 0, true, iterationnum);
	}

	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof IntensityFeature) ) return false;
    	IntensityFeature con = (IntensityFeature)o;
    	
    	return super.equals(con);
	}

	@Override
	public int hashCode() {
		return super.hashCode();
	}

	@Override
	public boolean isSatisfiedBy(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol,  Tolerance pmtol, int iterationNum) {
		return super.isSatisfiedBy(bp, iterationNum, spar, ppar);
	}

	@Override
	public String toFileString() {
		return "N\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString() +"\t" +this.getBasePeakCharge();
	}
	
	@Override
	public String toString() {
		return "Null - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge();
	}


}
