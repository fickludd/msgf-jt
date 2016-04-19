package uninovoOld.parameters;

import java.util.ArrayList;
import java.util.HashSet;

import msutil.Peak;
import msutil.Spectrum;

public class PeakParameter {
	private static int intensityRatioLevelNum = 5 * 2;
	static private int intensityLevelNum = 10;
	static private int partitionNum = 4; // 0 => partition num varies
	static private int ionCurrentRatioNum = 5;
	static private float[][] partitionMzs = null;
	static private Spectrum currentSpec = null;
	static private float ionCurrent = -1, currentHighestPeakIntensity = 0;
	
	private int basePeakGroup, basePeakPartition, basePeakPartition2=-1;
	private int ionCurrentRatio = -1, iterationNum = -1;
	
	private PeakParameter(int basePeakGroup, int basePeakPartition, int basePeakPartition2, int ionCurrentRatio){
		this(basePeakGroup, basePeakPartition, ionCurrentRatio);
		this.basePeakPartition2 = basePeakPartition2;
	}
	
	private PeakParameter(int basePeakGroup, int basePeakPartition, int ionCurrentRatio){
		this.basePeakGroup = basePeakGroup;
		this.basePeakPartition = basePeakPartition;
		this.ionCurrentRatio = ionCurrentRatio;
	}
	
	public PeakParameter(Peak bp, Spectrum spec, int iterationNum){
		if(currentSpec == null || currentSpec != spec){
			currentSpec = spec;
			ionCurrent = 0;
			for(Peak p : spec){
				ionCurrent += p.getIntensity();		
				if(p.getRank() == 1) currentHighestPeakIntensity = p.getIntensity();
				//	currentHighestPeakIntensity = Math.max(currentHighestPeakIntensity, p.getIntensity());
			}
		}
		
		this.iterationNum = iterationNum;
		
		if(intensityLevelNum == 1) basePeakGroup = 0;
		else if(iterationNum == 0){
			if(basePeakPartition == 1)
				basePeakGroup = Math.min(intensityLevelNum - 1, bp.getRank()/(100/(intensityLevelNum-1)));
			else{
				double x = Math.pow(100, 1.0/(intensityLevelNum - 1));			
				basePeakGroup = Math.min(intensityLevelNum - 1, (int)(Math.log(bp.getRank())/Math.log(x)));
			}//System.out.println(basePeakGroup + "\t" + bp.getRank());
		}else{
			//basePeakGroup = intensityLevelNum - 1 - Math.round((intensityLevelNum-2) * bp.getIntensity()/currentHighestPeakIntensity);
			//double x = Math.pow(100, 1.0/intensityLevelNum);
			//double y = Math.max(1-bp.getIntensity()/currentHighestPeakIntensity, 0.01);
			
			//basePeakGroup = Math.max(0, intensityLevelNum+(int)(Math.log(y)/Math.log(x)));
			//basePeakGroup = Math.min(basePeakGroup, intensityLevelNum - 1);
			//basePeakGroup = Math.min(intensityLevelNum - 1, (int)(Math.log(bp.getRank())/Math.log(x)));
			basePeakGroup = intensityLevelNum - 1 - Math.round((intensityLevelNum-2) * bp.getIntensity()/currentHighestPeakIntensity);
			
		}
		if(partitionNum <= 0){
			basePeakPartition = 1;
			for(int charge = Math.max(1, Math.min(3, spec.getCharge())); charge >=1; charge--){//TODO make '3' as a variable?
				if(bp.getMz() < maxMzWith(charge, spec)){
					basePeakPartition = charge;
					break;
				}
			}
		}else{
			float partitionSize = spec.getParentMass() / partitionNum;
			basePeakPartition =  (int)Math.min((bp.getMz() / partitionSize), partitionNum-1);
		}
			//if(bp.getMz() > spec.getParentMass()/4*3) basePeakPartition = spec.getCharge() + 1;
			//float partitionSize = spec.getParentMass() / partitionNum;
			//basePeakPartition =  (int)Math.min((bp.getMz() / partitionSize) + 1, partitionNum);
		//}
		
		if(partitionMzs != null && partitionMzs[spec.getCharge()] != null){
			int i=0;
			for(i=0; i< partitionMzs[spec.getCharge()].length;i++){
				if(bp.getMz() < partitionMzs[spec.getCharge()][i])
					break;
			}
			//System.out.println(spec.getPeptideMass() + "\t" + i);
			basePeakPartition2 = i;
		}
		
		if(ionCurrentRatioNum>1 && iterationNum == 0 && basePeakGroup == 0){
			ionCurrentRatio = Math.min(ionCurrentRatioNum-1, Math.round(4 * bp.getIntensity()/ionCurrent * (ionCurrentRatioNum-1)));			
		}else ionCurrentRatio = -1;
		
	}
	
	public String toFileString(){
		return basePeakGroup + " " + basePeakPartition + (basePeakPartition2>=0? " " + basePeakPartition2 : "") + " " + ionCurrentRatio;
	}
	
	public String toString(){
		return "Group: " + basePeakGroup + " Partition: " + basePeakPartition + " Partition2: " + basePeakPartition2 + " IonCurrentratio: "+ionCurrentRatio;
	}
	
	public int getBasePeakGroupNum() { return basePeakGroup; }
	public int getBasePeakPartitionNum() { return basePeakPartition; }
	public int getBasePeakPartitionNum2() { return basePeakPartition2; }
	public int getIonCurrentRatio() { return ionCurrentRatio; }
	
	public boolean equals(Object o){
		if(this == o) return true;
    	if ( !(o instanceof PeakParameter) ) return false;
    	PeakParameter ppar = (PeakParameter)o;
    	
    	return this.basePeakGroup == ppar.basePeakGroup 
    		&& this.basePeakPartition == ppar.basePeakPartition
    		&& this.basePeakPartition2 == ppar.basePeakPartition2
    		&& this.ionCurrentRatio == ppar.ionCurrentRatio;
	}
	
	public int hashCode(){
		return basePeakGroup + (basePeakPartition << 8) + (basePeakPartition2 << 4) + (ionCurrentRatio << 6);
	}
	
	static public void setGroupNum(int g) { intensityLevelNum = g; }
	
	static public void setPartitionNum(int p) {partitionNum = p;}
	
	static public void setIonCurrentRatioNum(int r){ ionCurrentRatioNum = r;}
	
	static public void setPartitionMzs(float[] partitionMzs, int charge){
		if(PeakParameter.partitionMzs == null) PeakParameter.partitionMzs = new float[100][];
		PeakParameter.partitionMzs[charge] = partitionMzs;
	//	System.out.println(writePartitionMzs(charge));
	}
	
	static public int getMaxGroupNum() {return intensityLevelNum;}
	
	static public PeakParameter parsePeakParameter(String s){
		String[] token = s.split(" ");
		if(token.length == 3)
			return new PeakParameter(Integer.parseInt(token[0]), Integer.parseInt(token[1]), Integer.parseInt(token[2]));
		else 
			return new PeakParameter(Integer.parseInt(token[0]), Integer.parseInt(token[1]), Integer.parseInt(token[2]), Integer.parseInt(token[3]));
	}
	
	static public float maxMzWith(int charge, Spectrum spec){		
		if(charge == 1) return Float.MAX_VALUE;
		
		return (float) (spec.getParentMass()/charge);
	}
	
	static public HashSet<Integer> getAllPartitionNum(int specCharge, int iterationNum){
		HashSet<Integer> bp = new HashSet<Integer>();
		for(PeakParameter pp : PeakParameter.getAllBasePeakParameters(specCharge, iterationNum)){
			bp.add(pp.getBasePeakPartitionNum());
		}
		return bp;
	}
	
	static public ArrayList<PeakParameter> getAllBasePeakParameters(int charge, int iterationNum){
		ArrayList<PeakParameter> pars = new ArrayList<PeakParameter>();
		int maxPartitionNum = partitionNum-1;
		int minPartitionNum = 0;
		if(partitionNum <= 0){
			maxPartitionNum = Math.min(3, charge);//TODO
			minPartitionNum = 1;
		}
		
		for(int p=minPartitionNum; p<=maxPartitionNum; p++){
			for(int g=0; g<intensityLevelNum; g++){
				if(partitionMzs != null && partitionMzs[charge] != null){
					for(int q=0; q<=partitionMzs[charge].length; q++){
						if(ionCurrentRatioNum>1 && iterationNum == 0 && g == 0){
							for(int i=0;i<ionCurrentRatioNum; i++)
								pars.add(new PeakParameter(g, p, q, i));
						}else pars.add(new PeakParameter(g, p, q, -1));
					}
				}
				else{
					if(ionCurrentRatioNum>1 && iterationNum == 0 && g == 0){
						for(int i=0;i<ionCurrentRatioNum; i++)
							pars.add(new PeakParameter(g, p, i));
					}else
						pars.add(new PeakParameter(g, p, -1));
				}
			}
		}
		return pars;
	}
	
	static public void setPeakIntensityRatioNum(int n){ intensityRatioLevelNum = n;}
	
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
	
	static public void refresh() {currentSpec = null;}
	
	static public ArrayList<Integer> getAllPeakIntensityRatioNums(){
		ArrayList<Integer> ratios = new ArrayList<Integer>();

		for(int r = -intensityRatioLevelNum/2; r<= intensityRatioLevelNum/2; r++)
			ratios.add(r);		
		return ratios;
	}

	
}
