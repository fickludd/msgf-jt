package uninovoOld.features;


import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Peak;
import msutil.Spectrum;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;
 
public class LinkingFeature extends Feature{
	private AminoAcidSet aaSet;
	
	//static private HashMap<Peak, HashMap<Integer, ArrayList<Integer>>> peakRatioMap = null;
	static private Spectrum currentSpec = null;
	static private int currentIterationNum = -1;
	
	public LinkingFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge, boolean isPresent,
			int iterationNum, int peakIntensityRatio, AminoAcidSet aaSet) {
		super(spar, ppar, basePeakCharge, peakIntensityRatio, isPresent, iterationNum);
		this.setPeakIntensityRatio(peakIntensityRatio);
		this.aaSet = aaSet;
	}

	private Peak currentPeak = null;
	private boolean currentHold = true;
	@Override
	public boolean holdsFor(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol, Tolerance pmtol, int iterationNum) {
		if(currentSpec!=null && spec.getParentMass() == currentSpec.getParentMass() && spec.equals(currentSpec) && bp == currentPeak){
			return currentHold;
		}
		
		boolean hold = false;

		if(super.holdsFor(bp, iterationNum, spar, ppar)){	
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
				//		ratios.add(ratio);
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
	public String toString() {
		return "Bridging - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge() + " Ratio: " + this.getPeakIntensityRatio() + " Present: " + this.isPresent();
	}

	@Override
	public String toFileString() {
		return "B\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString()+"\t" +this.getBasePeakCharge() + "\t" + this.getPeakIntensityRatio() + "\t" + (this.isPresent()? "p" : "a");
	}
	
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




}
