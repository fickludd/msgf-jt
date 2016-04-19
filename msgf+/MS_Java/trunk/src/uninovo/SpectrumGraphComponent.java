/***************************************************************************
  * Title:          
  * Author:        Kyowon Jeong
  * Last modified:  
  *
  * Copyright (c) 2008-2009 The Regents of the University of California
  * All Rights Reserved
  * See file LICENSE for details.
  ***************************************************************************/
package uninovo;

import java.io.IOException;
import java.util.Comparator;

import uninovo.parser.BufferedLineReader;
import uninovo.util.Enzyme;
import uninovo.util.Peptide;
import uninovo.util.Tolerance;

/**
 * The Class SpectrumGraphComponent is a super class for node and edge
 */
public abstract class SpectrumGraphComponent  implements Comparable<SpectrumGraphComponent>
{
	
	/**
	 * The Class LRComparator.
	 */
	public static class LRComparator implements Comparator<SpectrumGraphComponent>{ // 
		/**
		 * Gets the lr comparator
		 *
		 * @return the lR comparator
		 */
		static public LRComparator get(){
			return new LRComparator();
		}
		
		@Override
		public int compare(SpectrumGraphComponent arg0, SpectrumGraphComponent arg1) {
			return new Float(arg0.LR).compareTo(new Float(arg1.LR));
		}
	}
	
	/** The charge ranges per fragmentation method index.. */
	static protected int[][] chargeRanges;
	
	/** The edge accuracies. */
	static protected float[][][][][][] edgeAccuracies; // typeIndex, charge, minAANum, accL, accR, mass devi
	
	/** The edge deviation accuracies. */
	static protected float [][][][][] edgeDeviationAccuracies; // typeIndex, charge correct/incorrect, minAA, deviation index, 
	
	/** The enzyme. */
	static protected Enzyme enzyme = null;
	
	/** The ion offset. */
	static protected float[][][] ionOffset;
	
	/** The ion weights. */
	static protected float[][][][] ionWeights;
	
	/** The Constant maxMassGroupNum. */
	static protected final int maxMassGroupNum = 8;
	
	/** The mass group number. */
	static protected int massGroupNumber = maxMassGroupNum; 
	
	/** The Constant maxFragmentationMethodIndexNumber. */
	static final int maxFragmentationMethodIndexNumber = 10;
	
	
	/** The node no peak accuracies. */
	static  protected float[][][] nodeNoPeakAccuracies;
	
	/** The node null accuracies. */
	static protected float[][][] nodeNullAccuracies;
	
	/** The node posterior probabilities. */
	static protected float[][][][] nodePosteriorProbabilities; // typeIndex, charge c/ic minAANum
	
	/** The MS1 tols per fragmentation method index. */
	static protected Tolerance[] pmtols;
	
	/** The MS2 tols per fragmentation method index. */
	static protected Tolerance[] tols;
	
	/**
	 * Gets the bounded value (upper and lower bounded).
	 *
	 * @param obj the obj
	 * @param up the up
	 * @param down the down
	 * @return the bounded value
	 */
	protected static double getBoundedValue(double obj, double up, double down){
		return Math.max(Math.min(obj, up), down);
	}
	
	/**
	 * Read the statistics in the parameter file
	 *
	 * @param para the parameter file
	 * @param fragmentationMethodIndex the fragmentationMethodIndex
	 */
	static public void read(String para, int fragmentationMethodIndex){
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(para);
			String s;
			int mode = -100;
			int index = 0;
			int index2 = 0;
			int charge = 0;
			if(chargeRanges == null)chargeRanges = new int[100][2];
			chargeRanges[fragmentationMethodIndex][0] = 1000;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("##NODE&EDGEPARAMETERS##")) continue;
				
				if(s.startsWith("#PMTOL")){
					if(pmtols == null) pmtols = new Tolerance[100];
					pmtols[fragmentationMethodIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
					mode = -100;
					continue;
				}
				if(s.startsWith("#TOL")){
					if(tols == null) tols = new Tolerance[100];
					tols[fragmentationMethodIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
					mode = -100;
					continue;
				}
				
				if(s.startsWith("#IONWEIGHT")){
					charge = Integer.parseInt(s.split("\t")[1]);
					chargeRanges[fragmentationMethodIndex][0] = Math.min(chargeRanges[fragmentationMethodIndex][0] , charge);
					chargeRanges[fragmentationMethodIndex][1] = Math.max(chargeRanges[fragmentationMethodIndex][1] , charge);
					continue;
				}else if(s.startsWith("#NULLPROB")){
					if(nodeNullAccuracies == null){
						nodeNullAccuracies = new float[maxFragmentationMethodIndexNumber][100][];
				
					}
					String[] t = s.split("\t");
					charge = Integer.parseInt(t[1]);
					massGroupNumber = t.length-2;
					
					if(nodeNullAccuracies[fragmentationMethodIndex][charge] == null)
						nodeNullAccuracies[fragmentationMethodIndex][charge] = new float[massGroupNumber];
					
					
					for(int i=2;i<t.length;i++)
						nodeNullAccuracies[fragmentationMethodIndex][charge][i-2] = Float.parseFloat(t[i]);
					
					continue;
				}else if(s.startsWith("#NOPEAKPROB")){
					if(nodeNoPeakAccuracies == null){
						nodeNoPeakAccuracies = new float[maxFragmentationMethodIndexNumber][100][];
					}
					String[] t = s.split("\t");
					charge = Integer.parseInt(t[1]);
					if(nodeNoPeakAccuracies[fragmentationMethodIndex][charge] == null)
						nodeNoPeakAccuracies[fragmentationMethodIndex][charge] = new float[maxMassGroupNum];
					
					for(int i=2;i<t.length;i++)
						nodeNoPeakAccuracies[fragmentationMethodIndex][charge][i-2] = Float.parseFloat(t[i]);
					
					continue;
				}else if(s.startsWith("#OFF")){
					mode = 0;
					if(ionOffset == null){
						ionOffset = new float[maxFragmentationMethodIndexNumber][100][];
					}
					if(ionOffset[fragmentationMethodIndex][charge] == null)
						ionOffset[fragmentationMethodIndex][charge] = new float[Integer.parseInt(s.split("\t")[2])];
					
					index = Integer.parseInt(s.split("\t")[1]);
					continue;
				}else if(s.startsWith("#WEIGHTS")){
					mode = 1;
					if(ionWeights == null){
						ionWeights = (new float[maxFragmentationMethodIndexNumber][100][][]);
					}
					
					if(ionWeights[fragmentationMethodIndex][charge] == null)
						ionWeights[fragmentationMethodIndex][charge] = new float[Integer.parseInt(s.split("\t")[3])][];
	
					index = 0;
					index2 = Integer.parseInt(s.split("\t")[1]);
					ionWeights[fragmentationMethodIndex][charge][index2] = new float[Integer.parseInt(s.split("\t")[2])];
					continue;
				}else if(s.startsWith("#EDGEACCURACY")){
					mode = 3;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[5]);
					if(edgeAccuracies == null)
						edgeAccuracies = new float[maxFragmentationMethodIndexNumber][100][][][][];
					
					if(edgeAccuracies[fragmentationMethodIndex][charge] == null)
						edgeAccuracies[fragmentationMethodIndex][charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])][Integer.parseInt(s.split("\t")[4])];
					continue;
				}				
				else if(s.startsWith("#DEVIATIONPROB1")){
					String[] token = s.split("\t");
					
					if(edgeDeviationAccuracies == null)
						edgeDeviationAccuracies = new float[maxFragmentationMethodIndexNumber][100][2][][];					
					charge = Integer.parseInt(token[1]);
					
					if(edgeDeviationAccuracies[fragmentationMethodIndex][charge][0] == null){
						edgeDeviationAccuracies[fragmentationMethodIndex][charge][0] = new float[token.length-1][Edge.massDeviationIndexNum];
						edgeDeviationAccuracies[fragmentationMethodIndex][charge][1] = new float[token.length-1][Edge.massDeviationIndexNum];
					}
					
					for(int k=0;k<token.length-2;k++){
						String[] nums = token[k+2].split(":");
						for(int i=0;i<nums.length;i++){
							edgeDeviationAccuracies[fragmentationMethodIndex][charge][0][k+1][i] = Float.parseFloat(nums[i]);
						}
					}
					
					continue;
				}else if(s.startsWith("#DEVIATIONPROB2")){
					String[] token = s.split("\t");
					charge = Integer.parseInt(token[1]);
					for(int k=0;k<token.length-2;k++){
						String[] nums = token[k+2].split(":");
						for(int i=0;i<nums.length;i++){
							edgeDeviationAccuracies[fragmentationMethodIndex][charge][1][k+1][i] = Float.parseFloat(nums[i]);
						}
					}
				}else if(s.startsWith("#NODEPOSTERIOR")){
					if(nodePosteriorProbabilities == null){
						nodePosteriorProbabilities = new float[maxFragmentationMethodIndexNumber][100][2][];
					}
					String[] token = s.split("\t");
					
					for(int k=0;k<token.length-2;k++){
						String[] nums = token[k+2].split(":");
						nodePosteriorProbabilities[fragmentationMethodIndex][charge][k] = new float[nums.length+1];
						for(int i=0;i<nums.length;i++){
							nodePosteriorProbabilities[fragmentationMethodIndex][charge][k][i+1] = Float.parseFloat(nums[i]);
						}
					}
				}			
			else if(s.startsWith("#END")){
					mode = -100;
					continue;
				}
				
				if(mode == -1){
					continue;
				}else if(mode == 0){
					ionOffset[fragmentationMethodIndex][charge][index] = Float.parseFloat(s);
					continue;
				}else if(mode == 1){
					ionWeights[fragmentationMethodIndex][charge][index2][index] = Float.parseFloat(s);
					index++;
				}else if(mode == 3){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								for(int j=0; j<t.split(":").length; j++){
									String u = t.split(":")[j];
									if(!u.isEmpty())
										edgeAccuracies[fragmentationMethodIndex][charge][index][index2][index3][j] = Float.parseFloat(u);
								}
								index3++;
							}
						}
						index2++;
					}
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Sets the enzyme.
	 *
	 * @param enzyme the new enzyme
	 */
	public static void setEnzyme(Enzyme enzyme) {SpectrumGraphComponent.enzyme = enzyme;}
	
	/**
	 * Sets the mass group number
	 *
	 * @param s the mass group number
	 */
	static public void setMassGroupNumber(int s){ massGroupNumber = s; }
	
	/**
	 * Use binning or not.
	 *
	 * @param tol the MS2 tol
	 * @return true, if using binning
	 */
	static public boolean useBinning(Tolerance tol){
		if(tol.getToleranceAsDa(100) >= 0.01f) return true;
		return false;
	}
	
	/** The accuracy. */
	protected float accuracy;
	
	/** The charge. */
	protected int charge;
	
	/** The fragmentation Method Index. */
	protected int fragmentationMethodIndex;
	
	/** The lr score. */
	protected float LR = 0;
	
	/** The mass. */
	protected float mass;
	
	/** The no peak accuracy. */
	protected float noPeakAccuracy;
	
	/** The null accuracy (average accuracy)*/
	protected float nullAccuracy;
	
	/** The MS2 tol. */
	protected Tolerance tol;
	
	/**
	 * Instantiates a new spectrum graph component.
	 *
	 * @param mass the mass
	 * @param charge the charge
	 * @param fragmentationMethodIndex the fragmentationMethodIndex
	 */
	public SpectrumGraphComponent(float mass, int charge, int fragmentationMethodIndex){
		this.mass = mass;
		this.charge = charge;
		this.fragmentationMethodIndex = fragmentationMethodIndex;
	}
	
	@Override
	public int compareTo(SpectrumGraphComponent o) {
		return new Float(this.getMass()).compareTo(new Float(o.getMass()));
	}
	
	@Override
	abstract public boolean equals(Object obj);
	
	/**
	 * Gets the accuracy.
	 *
	 * @return the accuracy
	 */
	public float getAccuracy(){
		return accuracy;
	}
	
	/**
	 * Gets the lr.
	 *
	 * @return the lr
	 */
	abstract public float getLR();
	
	/**
	 * Gets the lR score.
	 *
	 * @return the integer lR score
	 */
	public int getLRScore(){
		return Math.round(getLR());
	}
	
	/**
	 * Gets the mass.
	 *
	 * @return the mass
	 */
	public float getMass(){
		return mass;
	}
	
	/**
	 * Gets the no peak accuracy.
	 *
	 * @return the no peak accuracy
	 */
	public float getNoPeakAccuracy(){
		return noPeakAccuracy;
	}
	
	/**
	 * Gets the null accuracy.
	 *
	 * @return the null accuracy
	 */
	public float getNullAccuracy(){
		return nullAccuracy;
	}
	
	@Override
	abstract public int hashCode();
	
	/**
	 * Checks if is correct.
	 *
	 * @param pep the pep
	 * @return true, if is correct
	 */
	public boolean isCorrect(Peptide pep){
		return isCorrect(pep, 0);
	}
	
	/**
	 * Checks if is correct.
	 *
	 * @param pep the peptide
	 * @param offset the offset
	 * @return true, if is correct
	 */
	abstract public boolean isCorrect(Peptide pep, float offset);
	
	/**
	 * Checks if is enzymatic.
	 *
	 * @return true, if is enzymatic
	 */
	abstract public boolean isEnzymatic();
	
	/**
	 * Checks if is within tolerance.
	 *
	 * @param mass the mass
	 * @param pm the pm
	 * @return true, if is within tolerance
	 */
	public boolean isWithinTolerance(float mass, float pm){
		return isWithinTolerance(mass, pm, null);
	}
	
	/**
	 * Checks if is within tolerance.
	 *
	 * @param mass the mass
	 * @param pm the pm
	 * @param tol the tol
	 * @return true, if is within tolerance
	 */
	public boolean isWithinTolerance(float mass, float pm, Tolerance tol){
		float diff = Math.abs(mass - getMass());
		if(tol == null) tol = this.tol;
		float t = tol.getToleranceAsDa(Math.max(pm - mass, mass));
		return diff <=t;
	} 
	
	/**
	 * Sets the accuracy.
	 *
	 * @param accuracy the new accuracy
	 */
	public void setAccuracy(float accuracy){
		this.accuracy = accuracy;
	}
	
	/**
	 * Sets the lr.
	 *
	 * @param LR the new lr
	 */
	public void setLR(float LR){
		this.LR = LR;
	}
	
	/**
	 * Sets the lr score and accuracy.
	 */
	protected abstract void setLRandAccuracy();
}
