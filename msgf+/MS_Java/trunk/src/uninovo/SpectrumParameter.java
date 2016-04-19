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

import uninovo.util.Spectrum;

/**
 * The Class SpectrumParameter  groups the spectra according to their mz and charge
 */
public class SpectrumParameter {
	
	/** The max spec mz range. */
	static private int[] maxSpecMzRange = null;
	
	/** The partition mzs. */
	static private float[][] partitionMzs = null;
	
	/** The spec mz range num. */
	static private int specMzRangeNum = 5;
	
	/**
	 * Gets the all spectrum parameters.
	 *
	 * @param charge the charge
	 * @return the all spectrum parameters
	 */
	static public ArrayList<SpectrumParameter> getAllSpectrumParameters(int charge){
		ArrayList<SpectrumParameter> pars = new ArrayList<SpectrumParameter>();
	
		int m = specMzRangeNum;
		if(maxSpecMzRange != null) m = Math.min(maxSpecMzRange[charge]+1, specMzRangeNum);
			
		for(int r=0; r<m; r++){
			pars.add(new SpectrumParameter(charge, r));
		}
		
		return pars;
	}
	
	/**
	 * Gets the max spec mz range.
	 *
	 * @param charge the charge
	 * @return the max spec mz range
	 */
	static public int getMaxSpecMzRange(int charge){
		if(maxSpecMzRange == null) return 0;
		return maxSpecMzRange[charge];
	}
	
	/**
	 * Gets the spec mz range.
	 *
	 * @param spec the spec
	 * @return the spec mz range
	 */
	static public int getSpecMzRange(Spectrum spec){ 
		if(specMzRangeNum == 1) return 0;
		if(partitionMzs!=null && partitionMzs[spec.getCharge()]!=null){
			int i=0;
			for(i=0; i< partitionMzs[spec.getCharge()].length;i++){
				if(spec.getPeptideMass() < partitionMzs[spec.getCharge()][i])
					break;
			}
			return i;
		}
		
		
			
		float len = spec.getPeptideMass() / 121.6f;
		if(specMzRangeNum == 2){
			if(len < 15) return 0;
			else return 1;
		}

		float m = (20f-9f)/(specMzRangeNum-2);
		
		for(int j=0; j<specMzRangeNum-1; j++){
			if(len <= 9f + j * m) return j;
		}
		return specMzRangeNum-1;
	}
	
	/**
	 * Gets the spec mz range num.
	 *
	 * @return the spec mz range num
	 */
	static public int getSpecMzRangeNum() {return specMzRangeNum;}
	
	/**
	 * Parses the spectrum parameter.
	 *
	 * @param s the s
	 * @return the spectrum parameter
	 */
	static public SpectrumParameter parseSpectrumParameter(String s){
		String[] token = s.split(" ");
		return new SpectrumParameter(Integer.parseInt(token[0]), Integer.parseInt(token[1]));
	}
	
	/**
	 * Sets the max spec mz range.
	 *
	 * @param charge the charge
	 * @param r the r
	 */
	static public void setMaxSpecMzRange(int charge, int r){
		if(maxSpecMzRange == null){
			maxSpecMzRange = new int[100];
		}
		maxSpecMzRange[charge] = r;
	}
	
	/**
	 * Sets the partition mzs.
	 *
	 * @param partitionMzs the partition mzs
	 * @param charge the charge
	 */
	static public void setPartitionMzs(float[] partitionMzs, int charge){
		if(SpectrumParameter.partitionMzs == null) SpectrumParameter.partitionMzs = new float[100][];
		SpectrumParameter.partitionMzs[charge] = partitionMzs;
	}
	
	/**
	 * Sets the spec mz range num.
	 *
	 * @param n the new spec mz range num
	 */
	static public void setSpecMzRangeNum(int n) {specMzRangeNum = n;}
	
	/**
	 * Write partition mzs.
	 *
	 * @param charge the charge
	 * @return the string
	 */
	static public String writePartitionMzs(int charge){
		if(partitionMzs == null || partitionMzs[charge] == null) return null;
		
		String out = null;

		boolean towrite = false;
		for(int i=0;i<partitionMzs[charge].length;i++){
			if(partitionMzs[charge][i]>0){
				towrite = true;
				break;
			}
		}
			
		if(towrite){
			out = "";
			for(int i=0;i<partitionMzs[charge].length;i++){
				if(partitionMzs[charge][i]>0)
					out += partitionMzs[charge][i]+"\t";
			}
		}
		return out;
		
	}
	
	/** The spec charge. */
	private int specCharge;
	
	/** The spec mz range. */
	private int specMzRange;
	
	/**
	 * Instantiates a new spectrum parameter.
	 *
	 * @param specCharge the spec charge
	 * @param specMzRange the spec mz range
	 */
	private SpectrumParameter(int specCharge, int specMzRange){
		this.specCharge = specCharge;
		this.specMzRange = specMzRange;
		if(maxSpecMzRange != null && maxSpecMzRange[specCharge]> 0)
			this.specMzRange = Math.min(maxSpecMzRange[specCharge], this.specMzRange);
	}
	
	/**
	 * Instantiates a new spectrum parameter.
	 *
	 * @param spec the spectrum
	 */
	public SpectrumParameter(Spectrum spec){
		this.specCharge = spec.getCharge();
		this.specMzRange = getSpecMzRange(spec);
		if(maxSpecMzRange != null && maxSpecMzRange[specCharge]> 0)
			this.specMzRange = Math.min(maxSpecMzRange[specCharge], this.specMzRange);
	}
	
	@Override
	public boolean equals(Object o){
		if(this == o) return true;
    	if ( !(o instanceof SpectrumParameter) ) return false;
    	SpectrumParameter spar = (SpectrumParameter)o;
    	return this.specCharge == spar.specCharge && this.specMzRange == spar.specMzRange;
    
	}
	
	/**
	 * Gets the spec charge.
	 *
	 * @return the spec charge
	 */
	public int getSpecCharge() { return specCharge; }
	
	/**
	 * Gets the spec mz range.
	 *
	 * @return the spec mz range
	 */
	public int getSpecMzRange() { return specMzRange; }
	
	@Override
	public int hashCode(){ return specCharge + (specMzRange << 8);}
	
	/**
	 * Checks if this parameter is for a given spectrum.
	 *
	 * @param spec the spectrum
	 * @return true, if is for
	 */
	public boolean isFor(Spectrum spec){
		return new SpectrumParameter(spec).equals(this);
	}
	
	/**
	 * To file string.
	 *
	 * @return the string
	 */
	public String toFileString() { return specCharge + " " + specMzRange;}
	
	@Override
	public String toString(){
		return "Charge: " + specCharge + " Mz: " + specMzRange;
	}
}
