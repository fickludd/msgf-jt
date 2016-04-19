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
import java.util.HashSet;

import uninovo.util.Peak;
import uninovo.util.Spectrum;

/**
 * The Class PeakParameter groups a pair of peaks according to their m/z, intensity, intensity ratio
 */
public class PeakParameter {
	
	/** The current spec. */
	static private Spectrum currentSpec = null;
	
	/** The intensity level total number. */
	static private int intensityLevelNum = 10;
	
	/** The intensity ratio level total number. */
	private static int intensityRatioLevelNum = 5 * 2;
	
	/** The current highest peak intensity. */
	static private float ionCurrent = -1, currentHighestPeakIntensity = 0;
	
	/** The ion current ratio total number. */
	static private int ionCurrentRatioNum = 5;
	
	/** The partition total number. */
	static private int partitionNum = 4; // 0 => partition num varies
	
	/**
	 * Gets the all base peak parameters.
	 *
	 * @param charge the charge
	 * @param iterationNum the iteration num
	 * @return the all base peak parameters
	 */
	static public ArrayList<PeakParameter> getAllBasePeakParameters(int charge, int iterationNum){
		ArrayList<PeakParameter> pars = new ArrayList<PeakParameter>();
		int maxPartitionNum = partitionNum-1;
		int minPartitionNum = 0;
		if(partitionNum <= 0){
			maxPartitionNum = Math.min(3, charge);
			minPartitionNum = 1;
		}
		
		for(int p=minPartitionNum; p<=maxPartitionNum; p++){
			for(int g=0; g<intensityLevelNum; g++){
				
				if(ionCurrentRatioNum>1 && iterationNum == 0 && g == 0){
					for(int i=0;i<ionCurrentRatioNum; i++)
						pars.add(new PeakParameter(g, p, i));
				}else
					pars.add(new PeakParameter(g, p, -1));
				
			}
		}
		return pars;
	}
	
	/**
	 * Gets the all partition num.
	 *
	 * @param specCharge the spec charge
	 * @param iterationNum the iteration num
	 * @return the all partition num
	 */
	static public HashSet<Integer> getAllPartitionNum(int specCharge, int iterationNum){
		HashSet<Integer> bp = new HashSet<Integer>();
		for(PeakParameter pp : PeakParameter.getAllBasePeakParameters(specCharge, iterationNum)){
			bp.add(pp.getBasePeakPartitionNum());
		}
		return bp;
	}
	
	/**
	 * Gets the all peak intensity ratio nums.
	 *
	 * @return the all peak intensity ratio nums
	 */
	static public ArrayList<Integer> getAllPeakIntensityRatioNums(){
		ArrayList<Integer> ratios = new ArrayList<Integer>();

		for(int r = -intensityRatioLevelNum/2; r<= intensityRatioLevelNum/2; r++)
			ratios.add(r);		
		return ratios;
	}
	
	/**
	 * Gets the max group num.
	 *
	 * @return the max group num
	 */
	static public int getMaxGroupNum() {return intensityLevelNum;}
	
	/**
	 * Gets the peak intensity ratio num.
	 *
	 * @param bp the bp
	 * @param cp the cp
	 * @param spec the spec
	 * @return the peak intensity ratio num
	 */
	static public int getPeakIntensityRatioNum(Peak bp, Peak cp, Spectrum spec){
		float bpi = bp.getIntensity();
		float cpi = cp.getIntensity();
		
		if(bpi == 0) bpi = Float.MIN_VALUE;
		if(cpi == 0) cpi = Float.MIN_VALUE;
		
		if(intensityRatioLevelNum == 1) return 0; 
		
		if(spec.getPrecursorPeak().equals(cp)) return -intensityRatioLevelNum/2;
		
		float x = bpi/cpi;
		
		if(x>=1){
			return  -(int) Math.floor((1-1/x)*(intensityRatioLevelNum/2));
		}
		else{
			return  (int) Math.floor((1-x)*(intensityRatioLevelNum/2)) + 1;
		}
	}
	
	/**
	 * Max mz with.
	 *
	 * @param charge the charge
	 * @param spec the spec
	 * @return the float
	 */
	static public float maxMzWith(int charge, Spectrum spec){		
		if(charge == 1) return Float.MAX_VALUE;
		
		return spec.getParentMass()/charge;
	}
	
	/**
	 * Parses the peak parameter.
	 *
	 * @param s the s
	 * @return the peak parameter
	 */
	static public PeakParameter parsePeakParameter(String s){
		String[] token = s.split(" ");
	
		return new PeakParameter(Integer.parseInt(token[0]), Integer.parseInt(token[1]), Integer.parseInt(token[2]));
		}
	
	/**
	 * Refresh.
	 */
	static public void refresh() {currentSpec = null;}
	
	/**
	 * Sets the group num.
	 *
	 * @param g the new group num
	 */
	static public void setGroupNum(int g) { intensityLevelNum = g; }
	
	/**
	 * Sets the ion current ratio num.
	 *
	 * @param r the new ion current ratio num
	 */
	static public void setIonCurrentRatioNum(int r){ ionCurrentRatioNum = r;}
	
	/**
	 * Sets the partition num.
	 *
	 * @param p the new partition num
	 */
	static public void setPartitionNum(int p) {partitionNum = p;}
	
	/**
	 * Sets the peak intensity ratio num.
	 *
	 * @param n the new peak intensity ratio num
	 */
	static public void setPeakIntensityRatioNum(int n){ intensityRatioLevelNum = n;}
	
	/** The base peak partition2. */
	private int basePeakGroup, basePeakPartition=-1;
	
	/** The ion current ratio. */
	private int ionCurrentRatio = -1;
	
	/**
	 * Instantiates a new peak parameter.
	 *
	 * @param basePeakGroup the base peak group
	 * @param basePeakPartition the base peak partition
	 * @param ionCurrentRatio the ion current ratio
	 */
	private PeakParameter(int basePeakGroup, int basePeakPartition, int ionCurrentRatio){
		this.basePeakGroup = basePeakGroup;
		this.basePeakPartition = basePeakPartition;
		this.ionCurrentRatio = ionCurrentRatio;
	}
	
	/**
	 * Instantiates a new peak parameter.
	 *
	 * @param bp the bp
	 * @param spec the spec
	 * @param iterationNum the iteration num
	 */
	public PeakParameter(Peak bp, Spectrum spec, int iterationNum){
		if(currentSpec == null || currentSpec != spec){
			currentSpec = spec;
			ionCurrent = 0;
			for(Peak p : spec){
				ionCurrent += p.getIntensity();		
				if(p.getRank() == 1) currentHighestPeakIntensity = p.getIntensity();
			}
		}	
		
		if(intensityLevelNum == 1) basePeakGroup = 0;
		else if(iterationNum == 0){
			if(basePeakPartition == 1)
				basePeakGroup = Math.min(intensityLevelNum - 1, bp.getRank()/(100/(intensityLevelNum-1)));
			else{
				double x = Math.pow(100, 1.0/(intensityLevelNum - 1));			
				basePeakGroup = Math.min(intensityLevelNum - 1, (int)(Math.log(bp.getRank())/Math.log(x)));
			}
		}else{
		basePeakGroup = intensityLevelNum - 1 - Math.round((intensityLevelNum-2) * bp.getIntensity()/currentHighestPeakIntensity);
			
		}
		if(partitionNum <= 0){
			basePeakPartition = 1;
			for(int charge = Math.max(1, Math.min(3, spec.getCharge())); charge >=1; charge--){
				if(bp.getMz() < maxMzWith(charge, spec)){
					basePeakPartition = charge;
					break;
				}
			}
		}else{
			float partitionSize = spec.getParentMass() / partitionNum;
			basePeakPartition =  (int)Math.min((bp.getMz() / partitionSize), partitionNum-1);
		}
		
		if(ionCurrentRatioNum>1 && iterationNum == 0 && basePeakGroup == 0){
			ionCurrentRatio = Math.min(ionCurrentRatioNum-1, Math.round(4 * bp.getIntensity()/ionCurrent * (ionCurrentRatioNum-1)));			
		}else ionCurrentRatio = -1;
		
	}
	
	@Override
	public boolean equals(Object o){
		if(this == o) return true;
    	if ( !(o instanceof PeakParameter) ) return false;
    	PeakParameter ppar = (PeakParameter)o;
    	
    	return this.basePeakGroup == ppar.basePeakGroup 
    		&& this.basePeakPartition == ppar.basePeakPartition
    		&& this.ionCurrentRatio == ppar.ionCurrentRatio;
	}
	
	/**
	 * Gets the base peak group num.
	 *
	 * @return the base peak group num
	 */
	public int getBasePeakGroupNum() { return basePeakGroup; }
	
	/**
	 * Gets the base peak partition num.
	 *
	 * @return the base peak partition num
	 */
	public int getBasePeakPartitionNum() { return basePeakPartition; }
	
	/**
	 * Gets the ion current ratio.
	 *
	 * @return the ion current ratio
	 */
	public int getIonCurrentRatio() { return ionCurrentRatio; }
	
	@Override
	public int hashCode(){
		return basePeakGroup + (basePeakPartition << 8) + (ionCurrentRatio << 6);
	}
	
	/**
	 * To file string.
	 *
	 * @return the string
	 */
	public String toFileString(){
		return basePeakGroup + " " + basePeakPartition + " " + ionCurrentRatio;
	}
	
	@Override
	public String toString(){
		return "Group: " + basePeakGroup + " Partition: " + basePeakPartition + " IonCurrentratio: "+ionCurrentRatio;
	}

	
}
