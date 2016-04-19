package uninovoOld.features;

import msgf.Tolerance;
import msutil.Peak;
import msutil.Spectrum;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;
import uninovoOld.train.InterPeakOffset;

public class OffsetFeature extends Feature{
	private InterPeakOffset gof;
	
	public OffsetFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge, boolean isPresent,
			int iterationNum, int peakIntensityRatio, InterPeakOffset gof) {
		super(spar, ppar, basePeakCharge, peakIntensityRatio, isPresent, iterationNum);
		assert(basePeakCharge == gof.getBaseCharge());
		this.setIOFF(gof);
	}


	private Peak currentPeak = null;
	private boolean currentHold = true;
	private Spectrum currentSpec = null;
	
	@Override
	public boolean holdsFor(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol, Tolerance pmtol, int iterationNum) {
		if(currentSpec!=null && spec.getParentMass() == currentSpec.getParentMass() && spec.equals(currentSpec) && bp == currentPeak){
			return currentHold;
		}
		
		boolean hold = false;

		if(super.holdsFor(bp, iterationNum, spar, ppar)){
			boolean match = false;
		
			for(Peak cp : getIOFF().getMatchingPeaks(bp, spec, tol, pmtol)){
				if(getPeakIntensityRatio() == PeakParameter.getPeakIntensityRatioNum(bp, cp, spec)){
					match = true;
					break;
				}
			}
			hold = match == this.isPresent();	
		}
		
		currentSpec = spec;
		currentHold = hold;
		currentPeak = bp;
		
		return hold;
	}

	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof OffsetFeature) ) return false;
    	OffsetFeature con = (OffsetFeature)o;
    	
    	return this.getPeakIntensityRatio() == con.getPeakIntensityRatio() && super.equals(con) && this.getIOFF().equals(con.getIOFF());
	}

	@Override
	public int hashCode() {
		return super.hashCode() * (this.getPeakIntensityRatio() + 48701) * this.getIOFF().hashCode();
	}

	@Override
	public String toString() {
		return "GOF - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge() + " Ratio: " + this.getPeakIntensityRatio() + " Present: " + this.isPresent() + " GOF: " + this.getIOFF();
	}
	
	@Override
	public String toFileString() {
		return "G\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString() +"\t" +this.getBasePeakCharge() + "\t" + this.getPeakIntensityRatio() + "\t" + (this.isPresent()? "p" : "a") + "\t" + this.getIOFF().toFileString();
	}
	
	public boolean isComplementary() {return getIOFF().isComplementary();}
	public int getChargeOffset() {return getIOFF().getChargeOffset();}
	
	static public OffsetFeature parseFileString(String s){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseSpectrumParameter(token[2]);
		PeakParameter ppar = PeakParameter.parsePeakParameter(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		int peakIntensityRatio = Integer.parseInt(token[5]);
		boolean isPresent = token[6].equals("p");
		InterPeakOffset gof = InterPeakOffset.parseFileString(token[7]);
		
		return new OffsetFeature(spar, ppar, basePeakCharge, isPresent, iterationNum, peakIntensityRatio, gof);
	}

	public void setIOFF(InterPeakOffset gof) {
		this.gof = gof;
	}

	public InterPeakOffset getIOFF() {
		return gof;
	}

	
}
