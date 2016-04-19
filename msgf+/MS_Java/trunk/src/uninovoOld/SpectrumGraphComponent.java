package uninovoOld;

import java.io.IOException;
import java.util.Comparator;

import msgf.Tolerance;
import msutil.Enzyme;
import msutil.Peptide;
import parser.BufferedLineReader;

public abstract class SpectrumGraphComponent  implements Comparable<SpectrumGraphComponent>
{
	protected float mass;
	protected float LR = 0;//-Float.MAX_VALUE;
	protected float accuracy;
	protected float nullAccuracy;
	protected float noPeakAccuracy;
	protected int charge;
	protected int typeIndex;
	protected Tolerance tol; 
	//public boolean isCorrect;
	
	static final int maxType = 10;
	static public final int maxSectionNum = 8;
	protected static int sectionNumber = maxSectionNum;
	static protected float[][][][] ionWeights;
	static protected float[][][] ionOffset;
	static protected Enzyme enzyme = null;
	static protected Tolerance[] tols;
	static protected Tolerance[] pmtols;
	static protected int[][] chargeRanges;
	
	static protected float[][][] nodeNullAccuracies;
	static  protected float[][][] nodeNoPeakAccuracies;
	
	static protected float[][][][][][] edgeAccuracies; // typeIndex, charge, minAANum, accL, accR, mass devi
	static protected float [][][][][] edgeDeviationAccuracies; // typeIndex, charge correct/incorrect, minAA, deviation index, 
	static protected float[][][][] nodePosteriorProbabilities; // typeIndex, charge c/ic minAANum
	
	static public boolean useBinning(Tolerance tol){
		if(tol.getToleranceAsDa(100) >= 0.01f) return true;
		return false;
	}
	
	public boolean isWithinTolerance(float mass, float pm){
		return isWithinTolerance(mass, pm, null);
	}
	
	public boolean isWithinTolerance(float mass, float pm, Tolerance tol){
		float diff = Math.abs(mass - getMass());
		if(tol == null) tol = this.tol;
		float t = tol.getToleranceAsDa(Math.max(pm - mass, mass));
		return diff <=t;
	}
	
	public SpectrumGraphComponent(float mass, int charge, int typeIndex){
		this.mass = mass;
		this.charge = charge;
		this.typeIndex = typeIndex;
	}
	
	protected abstract void setLRandAccuracy();
	
	protected static double getBoundedValue(double obj, double up, double down){
		return Math.max(Math.min(obj, up), down);
	}
	
	public int getLRScore(){
		return Math.round(getLR());
	}
	
	abstract public float getLR();
	
	public void setLR(float LR){
		this.LR = LR;
	}
	
	public float getAccuracy(){
		return accuracy;
	}
	
	public void setAccuracy(float accuracy){
		this.accuracy = accuracy;
	}
	
	public float getNullAccuracy(){
		return nullAccuracy;
	}
	
	public float getNoPeakAccuracy(){
		return noPeakAccuracy;
	}
	
	public float getMass(){
		return mass;
	}
	
	@Override
	public int compareTo(SpectrumGraphComponent o) {
		return new Float(this.getMass()).compareTo(new Float(o.getMass()));
	}
	
	@Override
	abstract public boolean equals(Object obj);
	
	abstract public boolean isCorrect(Peptide pep, float offset);
	
	public boolean isCorrect(Peptide pep){
		return isCorrect(pep, 0);
	}
	
	@Override
	abstract public int hashCode();
	
	abstract public boolean isEnzymatic();
	
	public static void setEnzyme(Enzyme enzyme) {SpectrumGraphComponent.enzyme = enzyme;} 
	
	static public void setSectionNumber(int s){ sectionNumber = s; }
	
	static public void read(String para, int typeIndex){
		//ArrayList<Integer> tm = new ArrayList<Integer>();
		//tm.add(specCharge); tm.add(typeIndex);
		//if(tmp.contains(tm)) return;
		//else tmp.add(tm);
		///if(edgeAccuracies!=null && edgeAccuracies[typeIndex]!=null && edgeAccuracies[typeIndex][specCharge] != null) return;
		
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(para);
			String s;
			int mode = -100;
			int index = 0;
			int index2 = 0;
			int charge = 0;
		//	int cur = 0;
			if(chargeRanges == null)chargeRanges = new int[100][2];
			chargeRanges[typeIndex][0] = 1000;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("##NODE&EDGEPARAMETERS##")) continue;
				
				if(s.startsWith("#PMTOL")){
					if(pmtols == null) pmtols = new Tolerance[100];
					pmtols[typeIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
					mode = -100;
					continue;
				}
				if(s.startsWith("#TOL")){
					if(tols == null) tols = new Tolerance[100];
					tols[typeIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
					mode = -100;
					continue;
				}
				
				if(s.startsWith("#IONWEIGHT")){
					charge = Integer.parseInt(s.split("\t")[1]);
					chargeRanges[typeIndex][0] = Math.min(chargeRanges[typeIndex][0] , charge);
					chargeRanges[typeIndex][1] = Math.max(chargeRanges[typeIndex][1] , charge);
					continue;
				}/*else if(s.startsWith("#ENZYMATICNODESACCURACY")){
					charge = Integer.parseInt(s.split("\t")[1]);
					if(accuracyForEnzymaticNodes == null){
						accuracyForEnzymaticNodes = new float[100][];
					}
					if(accuracyForEnzymaticNodes[charge] == null){
						accuracyForEnzymaticNodes[charge] = new float[Integer.parseInt(s.split("\t")[2])];
					}
					index = 0;
					mode = 4;
					continue;
				}*/else if(s.startsWith("#NULLPROB")){
					if(nodeNullAccuracies == null){
						nodeNullAccuracies = new float[maxType][100][];
					//	enzymaticNodesNullAccuracies = new float[100];
					}
					String[] t = s.split("\t");
					charge = Integer.parseInt(t[1]);
					sectionNumber = t.length-2;
					
					if(nodeNullAccuracies[typeIndex][charge] == null)
						nodeNullAccuracies[typeIndex][charge] = new float[sectionNumber];
					
					
					for(int i=2;i<t.length;i++)
						nodeNullAccuracies[typeIndex][charge][i-2] = Float.parseFloat(t[i]);
					
					continue;
				}else if(s.startsWith("#NOPEAKPROB")){
					if(nodeNoPeakAccuracies == null){
						nodeNoPeakAccuracies = new float[maxType][100][];
					}
					String[] t = s.split("\t");
					charge = Integer.parseInt(t[1]);
					if(nodeNoPeakAccuracies[typeIndex][charge] == null)
						nodeNoPeakAccuracies[typeIndex][charge] = new float[maxSectionNum];
					
					for(int i=2;i<t.length;i++)
						nodeNoPeakAccuracies[typeIndex][charge][i-2] = Float.parseFloat(t[i]);
					
					continue;
				}else if(s.startsWith("#OFF")){
					mode = 0;
					if(ionOffset == null){
						ionOffset = new float[maxType][100][];
					}
					if(ionOffset[typeIndex][charge] == null)
						ionOffset[typeIndex][charge] = new float[Integer.parseInt(s.split("\t")[2])];
					
					index = Integer.parseInt(s.split("\t")[1]);
					continue;
				}else if(s.startsWith("#WEIGHTS")){
					mode = 1;
					//	weightsForEachIon = new float[Integer.parseInt(s.split("\t")[1])][2];index = 0;
					if(ionWeights == null){
						ionWeights = (new float[maxType][100][][]);
					}
					
					if(ionWeights[typeIndex][charge] == null)
						ionWeights[typeIndex][charge] = new float[Integer.parseInt(s.split("\t")[3])][];
	
					index = 0;
					index2 = Integer.parseInt(s.split("\t")[1]);
					ionWeights[typeIndex][charge][index2] = new float[Integer.parseInt(s.split("\t")[2])];
					continue;
				}else if(s.startsWith("#EDGEACCURACY")){
					mode = 3;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[5]);
					if(edgeAccuracies == null)
						edgeAccuracies = new float[maxType][100][][][][];
					
					if(edgeAccuracies[typeIndex][charge] == null)
						edgeAccuracies[typeIndex][charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])][Integer.parseInt(s.split("\t")[4])];
					continue;
				}/*else if(s.startsWith("#IONPROBCOEFF")){
					mode = 5;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[6]);
					cur = Integer.parseInt(s.split("\t")[1]);
					if(ionProbCoeffs == null)
						ionProbCoeffs = new float[2][100][][][][];
					
					if(ionProbCoeffs[cur][charge] == null){
						ionProbCoeffs[cur][charge] = new float[Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])][Integer.parseInt(s.split("\t")[4])][Integer.parseInt(s.split("\t")[5])];	
					}
					continue;
				}else if(s.startsWith("#NULLEDGEACCURACY")){
					mode = 2;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[1]);
					
					if(edgeNullAccuracies == null)
						edgeNullAccuracies = new float[maxType][100];
					
					continue;
				}*/
				
				else if(s.startsWith("#DEVIATIONPROB1")){
					String[] token = s.split("\t");
					
					if(edgeDeviationAccuracies == null)
						edgeDeviationAccuracies = new float[maxType][100][2][][];
					//System.out.println(token.length);
					charge = Integer.parseInt(token[1]);
					
					if(edgeDeviationAccuracies[typeIndex][charge][0] == null){
						edgeDeviationAccuracies[typeIndex][charge][0] = new float[token.length-1][Edge.massDeviationIndexNum];
						edgeDeviationAccuracies[typeIndex][charge][1] = new float[token.length-1][Edge.massDeviationIndexNum];
					}
					
					for(int k=0;k<token.length-2;k++){
						String[] nums = token[k+2].split(":");
						for(int i=0;i<nums.length;i++){
							edgeDeviationAccuracies[typeIndex][charge][0][k+1][i] = Float.parseFloat(nums[i]);
						}
					}
					
					continue;
				}else if(s.startsWith("#DEVIATIONPROB2")){
					String[] token = s.split("\t");
					charge = Integer.parseInt(token[1]);
					for(int k=0;k<token.length-2;k++){
						String[] nums = token[k+2].split(":");
						for(int i=0;i<nums.length;i++){
							edgeDeviationAccuracies[typeIndex][charge][1][k+1][i] = Float.parseFloat(nums[i]);
						}
					}
				}else if(s.startsWith("#NODEPOSTERIOR")){
					if(nodePosteriorProbabilities == null){
						nodePosteriorProbabilities = new float[maxType][100][2][];
					}
					String[] token = s.split("\t");
					
					for(int k=0;k<token.length-2;k++){
						String[] nums = token[k+2].split(":");
						nodePosteriorProbabilities[typeIndex][charge][k] = new float[nums.length+1];
						for(int i=0;i<nums.length;i++){
							nodePosteriorProbabilities[typeIndex][charge][k][i+1] = Float.parseFloat(nums[i]);
						}
					}
				}
				
				
				/*else if(s.startsWith("#NODEACCURACY")){
				
					mode = 6;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[4]);
					if(nodeAccuracies == null)
						nodeAccuracies = new float[100][][][];
					
					if(nodeAccuracies[charge] == null)
						nodeAccuracies[charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])];
					continue;
				}*//*else if(s.startsWith("#ENZNODEACCURACY")){
					mode = 7;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[4]);
					if(enzymaticNodeAccuracies == null)
						enzymaticNodeAccuracies = new float[100][][][];
					
					if(enzymaticNodeAccuracies[charge] == null)
						enzymaticNodeAccuracies[charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])];
					continue;
				}else if(s.startsWith("#IONPROBPROB")){
					mode = 8;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[5]);
					if(ionProbProbabilities == null)
						ionProbProbabilities = new float[100][][][][];
					
					if(ionProbProbabilities[charge] == null)
						ionProbProbabilities[charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])][Integer.parseInt(s.split("\t")[4])];
					continue;
				}else if(s.startsWith("#ENZIONPROBPROB")){
					mode = 9;
					index = 0; index2 = 0;
					charge = Integer.parseInt(s.split("\t")[5]);
					if(enzymaticIonProbProbabilities == null)
						enzymaticIonProbProbabilities = new float[100][][][][];
					
					if(enzymaticIonProbProbabilities[charge] == null)
						enzymaticIonProbProbabilities[charge] = new float[Integer.parseInt(s.split("\t")[1])][Integer.parseInt(s.split("\t")[2])][Integer.parseInt(s.split("\t")[3])][Integer.parseInt(s.split("\t")[4])];
					continue;
				}*/else if(s.startsWith("#END")){
					mode = -100;
					continue;
				}
				
				if(mode == -1){
					continue;
				}else if(mode == 0){
					ionOffset[typeIndex][charge][index] = Float.parseFloat(s);
					continue;
				}else if(mode == 1){
					ionWeights[typeIndex][charge][index2][index] = Float.parseFloat(s);
					index++;
				}/*else if(mode == 2){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						edgeNullAccuracies[typeIndex][charge] = Float.parseFloat(s);
						
					}
				}*/else if(mode == 3){
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
										edgeAccuracies[typeIndex][charge][index][index2][index3][j] = Float.parseFloat(u);
								}
								index3++;
							}
						}
						index2++;
					}
				}/*else if(mode == 4){
					accuracyForEnzymaticNodes[charge][index] = Float.parseFloat(s);
					index++;						
				}*//*else if(mode == 5){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								String[] u = t.split(":");
								for(int l=0;l<u.length;l++){
									ionProbCoeffs[cur][charge][index][index2][index3][l] = Float.parseFloat(u[l]);
								}
								index3++;
							}
						}
						index2++;
					}
				}else if(mode == 6){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								nodeAccuracies[charge][index][index2][index3] = Float.parseFloat(t);
								index3++;
							}
						}
						index2++;
					}
				}*//*else if(mode == 7){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								enzymaticNodeAccuracies[charge][index][index2][index3] = Float.parseFloat(t);
								index3++;
							}
						}
						index2++;
					}
				}else if(mode == 8){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								String[] u = t.split(":");
								for(int l=0;l<u.length;l++){
									ionProbProbabilities[charge][index][index2][index3][l] = Float.parseFloat(u[l]);
								}
								index3++;
							}
						}
						index2++;
					}
				}else if(mode == 9){
					if(s.startsWith("##NUM")){
						index = Integer.parseInt(s.split("\t")[1]);
						index2 = 0;
						continue;
					}else{
						int index3=0;
					
						for(String t:s.split(" ")){
							if(!t.isEmpty()){
								String[] u = t.split(":");
								for(int l=0;l<u.length;l++){
									enzymaticIonProbProbabilities[charge][index][index2][index3][l] = Float.parseFloat(u[l]);
								}
								index3++;
							}
						}
						index2++;
					}
				}*/
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static class LRComparator implements Comparator<SpectrumGraphComponent>{ // 
		@Override
		public int compare(SpectrumGraphComponent arg0, SpectrumGraphComponent arg1) {
			return new Float(arg0.LR).compareTo(new Float(arg1.LR));
		}
		
		static public LRComparator get(){
			return new LRComparator();
		}
	}
}
