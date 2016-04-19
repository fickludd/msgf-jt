package uninovoOld.train;

import java.util.ArrayList;

import msgf.Tolerance;
import msutil.Constants;
import msutil.Peak;
import msutil.Spectrum;

public class InterPeakOffset {
	private float offset;
	private boolean isComplementary;
	private int chargeOffset;
	private int baseCharge;
	
	private InterPeakOffset(float offset, boolean isComplementary, int chargeOffset, int baseCharge){
		this.offset = offset;
		this.isComplementary = isComplementary;
		this.chargeOffset = chargeOffset;
		this.baseCharge = baseCharge;
	}
	
	public InterPeakOffset(float offset, boolean isComplementary, int chargeOffset, int baseCharge, float resolution){
		if(resolution > 0){
			this.offset = (Math.round(offset*Constants.INTEGER_MASS_SCALER_HIGH_PRECISION))/Constants.INTEGER_MASS_SCALER_HIGH_PRECISION;
		//	this.offset = Math.round(offset*Constants.INTEGER_MASS_SCALER*20) / Constants.INTEGER_MASS_SCALER / 20f;
			
		}
		else this.offset = offset;

		if(resolution >= 0.5) 
			this.offset = Math.round(offset*Constants.INTEGER_MASS_SCALER) / Constants.INTEGER_MASS_SCALER;
		
		this.isComplementary = isComplementary;
		this.chargeOffset = chargeOffset;
		this.baseCharge = baseCharge;
	}
	
	public float getOffset(){ return offset; }
	public boolean isComplementary() { return isComplementary; }
	public int getChargeOffset() { return chargeOffset; }
	public int getBaseCharge() { return baseCharge; }
	
	static public InterPeakOffset parseFileString(String s){ 
		String[] token = s.split(" ");
		
		return new InterPeakOffset(Float.parseFloat(token[0]), Integer.parseInt(token[1]) == 1,
				Integer.parseInt(token[2]), Integer.parseInt(token[3]));
	}
	
	public String toFileString(){
		String s = this.offset + " " + (this.isComplementary? 1 : 0) + " " + this.chargeOffset + " " + this.baseCharge;
		return s;
	}
	
	public String toString(){
		return "Off: " + this.offset + " Comp: " + this.isComplementary + " Charge: " + this.baseCharge + " ChargeOff: " + this.chargeOffset;
	}
	
	public int hashCode(){
		return this.baseCharge * (this.chargeOffset + 10) * (int)(this.offset) * (isComplementary? 1 : -1); 
	}
	
	public boolean equals(Object o){
    	if(this == o) return true;
    	if ( !(o instanceof InterPeakOffset) ) return false;
    	InterPeakOffset io = (InterPeakOffset)o;
    	return io.offset == this.offset && io.baseCharge == this.baseCharge
    		&& io.isComplementary == this.isComplementary && io.chargeOffset == this.chargeOffset;
    }
	
	public ArrayList<Peak> getMatchingPeaks(Peak bp, Spectrum spec, Tolerance tol, Tolerance pmtol){
		if(bp.getCharge() != baseCharge) return new ArrayList<Peak>();
		
		Tolerance to = tol;
		
		if(isComplementary){
			if(tol.getToleranceAsDa(100) < pmtol.getToleranceAsDa(100)) to = pmtol;
		}
		
		float t = to.getToleranceAsDa(bp.getMz() * bp.getCharge());
		
		//t /=  bp.getCharge();
		
		if(chargeOffset != 0)
			bp = PeakGenerator.getChargeChangedBasePeak(bp, bp.getCharge(), chargeOffset);
		if(isComplementary)
			bp = PeakGenerator.getComplementaryBasePeak(bp, bp.getCharge(), spec);
			
		float mz = bp.getMz() + offset / bp.getCharge();
		
		float s = to.getToleranceAsDa(mz * bp.getCharge());
		
		//s /=  bp.getCharge();
		
		if(isComplementary && t < 0.25f){
			float  u = pmtol.getToleranceAsDa(mz * bp.getCharge());
			//u /=  bp.getCharge();
			s = s>u? s : u;
			t += s;			
		}else	
			t = t>s? t : s;
		
		t = Math.min(t, 0.5f);
		ArrayList<Peak> ret = spec.getPeakListByMassRange(mz-t, mz + t);
		
		if(bp.getCharge() == spec.getCharge() && Math.abs(spec.getPrecursorPeak().getMz() - mz) < pmtol.getToleranceAsDa(mz))
			ret.add(spec.getPrecursorPeak());
	
		return ret;
	}
	
}
