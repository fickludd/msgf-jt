package uninovoOld.features;



import msgf.Tolerance;
import msutil.Peak;
import msutil.Spectrum;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;

public class IntensityFeature extends Feature{
	
	
	public IntensityFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge,
			int iterationnum) {
		super(spar, ppar, basePeakCharge, 0, true, iterationnum);
	}

	@Override
	public boolean holdsFor(Peak bp, Spectrum spec, SpectrumParameter spar, PeakParameter ppar, Tolerance tol,  Tolerance pmtol, int iterationNum) {
		return super.holdsFor(bp, iterationNum, spar, ppar);
	}

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
	public String toString() {
		return "Null - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge();
	}

	@Override
	public String toFileString() {
		return "N\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString() +"\t" +this.getBasePeakCharge();
	}
	
	static public IntensityFeature parseFileString(String s){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseSpectrumParameter(token[2]);
		PeakParameter ppar = PeakParameter.parsePeakParameter(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		
		return new IntensityFeature(spar,	ppar, basePeakCharge, iterationNum);
	}


}
