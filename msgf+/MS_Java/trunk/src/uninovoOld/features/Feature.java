package uninovoOld.features;

import java.util.ArrayList;
import java.util.HashMap;

import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.Composition;
import msutil.IonType;
import msutil.Peak;
import msutil.Spectrum;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;

public abstract class Feature implements Comparable<Feature>{
	private SpectrumParameter spar;
	private PeakParameter ppar;
	private int basePeakCharge;
	private int peakIntensityRatio;
	private boolean isPresent;
	private int iterationNum;
	private HashMap<IonType, Float> distribution = null;
	private HashMap<IonType, Float> number = null;// only for training
	
	private Feature intensityFeature = null;
	//TODO distribution shouldnt include precuror should include noisy 
	protected Feature(SpectrumParameter spar, PeakParameter ppar, int basePeakCharge, int peakIntensityRatio, boolean isPresent, int iterationnum){
		this.spar = spar;
		this.ppar = ppar;
		this.basePeakCharge = basePeakCharge;
		this.peakIntensityRatio = peakIntensityRatio;
		this.isPresent = isPresent;
		this.iterationNum = iterationnum;
		if(this.distribution == null) this.distribution = new HashMap<IonType, Float>();
	}

	public boolean isPresent() {return isPresent;}
	public int getIterationNum() {return iterationNum;}
	public int getBasePeakCharge() {return basePeakCharge;}

	public float getKLDivergenceFromIntensityFeature(){
		if(intensityFeature == null || this instanceof IntensityFeature) return 0;
		
		float kl = 0;
		
		for(IonType ion:this.distribution.keySet()){
			if(ion instanceof IonType.PrecursorIon) continue;
			float p1 = getIonProbability(ion);
			if(p1 <= 0) continue;
			float p2 = intensityFeature.getIonProbability(ion);	
			
			assert(p2 !=0);
			
			kl+= (float)(p1 * Math.log(p1/p2) / Math.log(2) );
		}

		return kl;
	}
	
	public SpectrumParameter getSpectrumParameter() {return spar;}
	public PeakParameter getBasePeakParameter() {return ppar;}
	
	public void registerNullCondition(Feature con){intensityFeature = con;}
	
	public Feature getIntensityFeature() {return intensityFeature;}
	
	public float getIonProbability(IonType ion){
	//	assert(distribution.isEmpty() || distribution.containsKey(IonType.NOISE));
		//System.out.println(distribution);
		if(this instanceof IntensityFeature){
			if(distribution.containsKey(ion))
				return distribution.get(ion);
			else return 0f;
		}else{
			float sum = 0;
			float nom = 0;
			for(IonType i : distribution.keySet()){
				float t = getProbability(i) * this.intensityFeature.getIonProbability(i);
				if(i.equals(ion)) nom = t;
				sum += t;
			}
			return nom/sum;
		}
	}
	
	public float getProbability(IonType ion){
	//	assert(distribution.isEmpty() || distribution.containsKey(IonType.NOISE));
		if(!(this instanceof IntensityFeature)){
			if(distribution.containsKey(ion))
				return distribution.get(ion);
			else return 0f;
		}
		return 0;
	}
	
	public int compareTo(Feature o) {
		return new Float(this.getKLDivergenceFromIntensityFeature()).compareTo(new Float(o.getKLDivergenceFromIntensityFeature()));
	}
	
	// ion should include Noise
	public void addIonCount(IonType ion){ 
		if(number == null) number = new HashMap<IonType, Float>();
		Float num = number.get(ion);
		if(num == null) num = 0f;
		num++;
		number.put(ion, num);
	}
	
	//returns sum
	public float calculateIonProbs(){
		if(number == null) return 0;
	//	if(!(this instanceof IntensityFeature)) return ;
		float sum = 0;
		for(IonType ion : number.keySet()){
			sum += number.get(ion);
		}
		
		for(IonType ion : number.keySet()){
			float p = 1;
			if(this instanceof IntensityFeature){
				p = number.get(ion)/sum;
			}else{
				p = number.get(ion)/this.intensityFeature.number.get(ion);
			}
			if(!(ion instanceof IonType.PrecursorIon)) distribution.put(ion, p);
		}
		return sum;
	}
	
	//noise!!
	public void setIonProbMap(HashMap<IonType, Float> ionProbMap){
		assert(ionProbMap.containsKey(IonType.NOISE));
		this.distribution = ionProbMap;
	}
	
	protected boolean equals(Feature con) {
		if(this == con) return true;
    	
    	return con.isPresent == this.isPresent &&
    		con.iterationNum == this.iterationNum &&
	    	con.basePeakCharge == this.basePeakCharge &&
	    	con.spar.equals(this.spar) &&
	    	con.ppar.equals(this.ppar);
	}
	
	protected boolean holdsFor(Peak bp, int iterationNum, SpectrumParameter spar, PeakParameter ppar){
		
		bp.setCharge(basePeakCharge);

		return iterationNum == this.iterationNum &&
		spar.equals(this.spar) &&
		ppar.equals(this.ppar);
	}
	
	public abstract boolean equals(Object o);
	
	public int hashCode(){
		int t = this.iterationNum + this.basePeakCharge << 3 + spar.hashCode() + ppar.hashCode();
		return this.isPresent? t : -t;
	}
	
	public abstract String toString();
	
	public abstract String toFileString();
	
	public abstract boolean holdsFor(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol, Tolerance pmtol, int iterationNum);

	static public float maxMzWithIon(IonType ion, Spectrum spec){
		if(ion.equals(IonType.NOISE)) return Float.MAX_VALUE;
		return ion.getMz((float) (spec.getParentMass() - Composition.H2O - AminoAcid.getStandardAminoAcid('G').getMass()));
	}
	
	static public float minMzWithIon(IonType ion){
		if(ion.equals(IonType.NOISE)) return -Float.MAX_VALUE;
		return ion.getMz(AminoAcid.getStandardAminoAcid('G').getMass());
	}
	
	static private float getUnNormalizedPeakIonProbability(IonType ion, ArrayList<Feature> cons, Peak bp, Spectrum spec){
		float p = 1;
		float py = 0;
		int n = 0;
		IntensityFeature nf = (IntensityFeature) cons.get(0).intensityFeature;
		
		py = nf.getIonProbability(ion);

		if(cons.size() == 1) return py;
		
		for(Feature con : cons){
			if(con instanceof IntensityFeature) continue;	
			float mult = 1;
			if(bp.getMz() < minMzWithIon(ion) || bp.getMz() > maxMzWithIon(ion, spec)) mult = 0;
			else mult = con.getProbability(ion); //TODO IonProb
			p*=mult;
			n ++;
		}
		
		if(py == 0) return 0;
		
		return (float) (p/Math.pow(py, n-1));
	}
	
	static private float getUnNormalizedPeakIonProbability2(IonType ion, ArrayList<Feature> cons, Peak bp, Spectrum spec){
		float p = 1;
		IntensityFeature nf = (IntensityFeature) cons.get(0).intensityFeature;
		
		p = nf.getIonProbability(ion);
		
		for(Feature con : cons){
			if(con instanceof IntensityFeature) continue;	
			float mult = 1;
			if(bp.getMz() < minMzWithIon(ion) || bp.getMz() > maxMzWithIon(ion, spec)) mult = 0;
			else mult = con.getProbability(ion);
			
			p*=mult;
		}

		return p;
	}

	static public HashMap<IonType, Float> getPeakIonProbabilityDistribution(Peak bp, ArrayList<Feature> cons, ArrayList<IonType> ions, HashMap<IonType, Float> rankProb, Spectrum spec){
		HashMap<IonType, Float> ret = new HashMap<IonType, Float>();
		 
		for(IonType ion : ions){ 
			ret.put(ion, getUnNormalizedPeakIonProbability2(ion, cons, bp , spec));
		}
		
		if(!ions.contains(IonType.NOISE)) ret.put(IonType.NOISE, getUnNormalizedPeakIonProbability2(IonType.NOISE, cons,  bp , spec));
		
		float sum = 0;
		for(IonType i : ret.keySet()){
			sum+=ret.get(i);
		}
	
		if(sum == 0){
			sum = 1;
			if(ret.containsKey(IonType.NOISE)) ret.put(IonType.NOISE, 1f);
		}
		for(IonType i : ret.keySet()){
			ret.put(i, ret.get(i)/sum);
		}

		
		return ret;
	}
	
	public void setPeakIntensityRatio(int peakIntensityRatio) {
		this.peakIntensityRatio = peakIntensityRatio;
	}

	public int getPeakIntensityRatio() {
		return peakIntensityRatio;
	}
	
}


