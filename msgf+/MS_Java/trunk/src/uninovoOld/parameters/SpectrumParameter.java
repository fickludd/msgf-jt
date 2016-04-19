package uninovoOld.parameters;

import java.util.ArrayList;

import msutil.Spectrum;

public class SpectrumParameter {
	private int specCharge;
	private int specMzRange;
	static private int specMzRangeNum = 5;
	static private float[][] partitionMzs = null;
	static private int[] maxSpecMzRange = null;
	
	private SpectrumParameter(int specCharge, int specMzRange){
		this.specCharge = specCharge;
		this.specMzRange = specMzRange;
		if(maxSpecMzRange != null && maxSpecMzRange[specCharge]> 0)
			this.specMzRange = Math.min(maxSpecMzRange[specCharge], this.specMzRange);
	}
	
	public SpectrumParameter(Spectrum spec){
		this.specCharge = spec.getCharge();
		this.specMzRange = getSpecMzRange(spec);
		if(maxSpecMzRange != null && maxSpecMzRange[specCharge]> 0)
			this.specMzRange = Math.min(maxSpecMzRange[specCharge], this.specMzRange);
	}
	
	public boolean isFor(Spectrum spec){
		return new SpectrumParameter(spec).equals(this);
	}
	
	static public void setMaxSpecMzRange(int charge, int r){
		if(maxSpecMzRange == null){
			maxSpecMzRange = new int[100];
		}
		maxSpecMzRange[charge] = r;
	}
	
	static public int getMaxSpecMzRange(int charge){
		if(maxSpecMzRange == null) return 0;
		return maxSpecMzRange[charge];
	}
	
	public String toFileString() { return specCharge + " " + specMzRange;}
	
	public String toString(){
		return "Charge: " + specCharge + " Mz: " + specMzRange;
	}
	
	public boolean equals(Object o){
		if(this == o) return true;
    	if ( !(o instanceof SpectrumParameter) ) return false;
    	SpectrumParameter spar = (SpectrumParameter)o;
    	return this.specCharge == spar.specCharge && this.specMzRange == spar.specMzRange;
    
	}
	
	public int hashCode(){ return specCharge + (specMzRange << 8);}
	
	public int getSpecCharge() { return specCharge; }
	public int getSpecMzRange() { return specMzRange; }
	
	static public void setPartitionMzs(float[] partitionMzs, int charge){
		if(SpectrumParameter.partitionMzs == null) SpectrumParameter.partitionMzs = new float[100][];
		SpectrumParameter.partitionMzs[charge] = partitionMzs;
	//	System.out.println(writePartitionMzs(charge));
	}
	
	static public int getSpecMzRange(Spectrum spec){ 
		if(specMzRangeNum == 1) return 0;
		if(partitionMzs!=null && partitionMzs[spec.getCharge()]!=null){
			int i=0;
			for(i=0; i< partitionMzs[spec.getCharge()].length;i++){
				if(spec.getPeptideMass() < partitionMzs[spec.getCharge()][i])
					break;
			}
			//System.out.println(spec.getPeptideMass() + "\t" + i);
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
	
	static public void setSpecMzRangeNum(int n) {specMzRangeNum = n;}
	static public int getSpecMzRangeNum() {return specMzRangeNum;}
	
	static public SpectrumParameter parseSpectrumParameter(String s){
		String[] token = s.split(" ");
		return new SpectrumParameter(Integer.parseInt(token[0]), Integer.parseInt(token[1]));
	}
	
	static public ArrayList<SpectrumParameter> getAllSpectrumParameters(int charge){
		ArrayList<SpectrumParameter> pars = new ArrayList<SpectrumParameter>();
	
		int m = specMzRangeNum;
		if(maxSpecMzRange != null) m = Math.min(maxSpecMzRange[charge]+1, specMzRangeNum);
			
		for(int r=0; r<m; r++){
			pars.add(new SpectrumParameter(charge, r));
		}
		
		return pars;
	}
}
