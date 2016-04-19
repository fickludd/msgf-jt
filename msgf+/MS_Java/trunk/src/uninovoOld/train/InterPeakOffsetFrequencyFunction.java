package uninovoOld.train;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import msgf.Tolerance;
import msutil.Peak;

public class InterPeakOffsetFrequencyFunction {
	static class InterPeakOffsetPeak  implements Comparable<InterPeakOffsetPeak>{
    	private InterPeakOffset gof;
    	private float y;
    	
    	InterPeakOffsetPeak(InterPeakOffset gof, float y){
    		this.gof = gof; this.y = y;
    	}
    	
    	public int compareTo(InterPeakOffsetPeak o) {
    		return new Float(this.y).compareTo(new Float(o.y));
    	}
    	
    	public float getProbability() { return y; }
    	public InterPeakOffset getInterPeakOffset() { return gof; }
    }
	
	static public final float MAX = 38;
	static public final float MIN = -38;
	
/*	static public ArrayList<NewGeneralizedOffsetPeak> getOffSetFrequencyFunction(HashMap<InterPeakOffset, Integer> offsetnums, float normalizer){
		return getOffSetFrequencyFunction(offsetnums, normalizer, 0);
	}
	
	static public ArrayList<NewGeneralizedOffsetPeak> getOffSetFrequencyFunction(HashMap<InterPeakOffset, Integer> offsetnums, float normalizer, float threshold){
		return getOffSetFrequencyFunction(offsetnums, normalizer, threshold, null);
	}
	*/
	static public ArrayList<InterPeakOffsetPeak> getOffSetFrequencyFunction(HashMap<InterPeakOffset, Integer> offsetnums, float normalizer, float threshold, String filename){
		ArrayList<InterPeakOffsetPeak> offsetPeaks = new  ArrayList<InterPeakOffsetPeak>();
		ArrayList<InterPeakOffsetPeak> offsetPeaksforOutput = new  ArrayList<InterPeakOffsetPeak>();
		if(offsetnums == null || offsetnums.isEmpty()) return offsetPeaks;
		
		if(threshold > 0 && threshold < 0.15f){
		
			ArrayList<Float> v = new ArrayList<Float>();
			float sum = 0, num = 0;
			for(InterPeakOffset gof : offsetnums.keySet()){
				v.add((float)offsetnums.get(gof)/normalizer);			
				sum+= (float)offsetnums.get(gof)/normalizer;
				num++;
			}
			Collections.sort(v);
		
			threshold = Math.min(threshold, v.get(v.size()/2)*7);
		}
		
		ArrayList<InterPeakOffsetPeak> offsetPeakstmp = new ArrayList<InterPeakOffsetPeak>();
		
		for(InterPeakOffset gof : offsetnums.keySet()){
			float prob = (float)offsetnums.get(gof)/normalizer;
			
			if(prob >= threshold)
				offsetPeakstmp.add(new InterPeakOffsetPeak(gof, prob));
			
			if(filename != null) offsetPeaksforOutput.add(new InterPeakOffsetPeak(gof, prob));
		}
		
		Collections.sort(offsetPeakstmp, Collections.reverseOrder());
		
		for(int i=0; i<Math.min(offsetPeakstmp.size(), 100); i++){
			offsetPeaks.add(offsetPeakstmp.get(i));
		}
		
		if(filename != null) Collections.sort(offsetPeaksforOutput);
		
		if(filename != null && !offsetPeaksforOutput.isEmpty()){
			try {
				HashSet<Integer> cs = new HashSet<Integer>();
				HashSet<Integer> cos = new HashSet<Integer>();
				HashSet<Boolean> comps = new HashSet<Boolean>();
				for(InterPeakOffsetPeak p : offsetPeaksforOutput){
					cs.add(p.gof.getBaseCharge());
					cos.add(p.gof.getChargeOffset());
					comps.add(p.gof.isComplementary());
				}
				
				for(int c : cs){
					for(int co : cos){
						for(boolean com : comps){
							boolean towrite = false;
							for(InterPeakOffsetPeak p : offsetPeaksforOutput){
								if(p.gof.getBaseCharge() == c && p.gof.getChargeOffset() == co && p.gof.isComplementary() == com)
									if(p.y >= threshold){
										towrite = true;
										break;
									}
							}
							
							if(!towrite) continue;
							
							PrintStream out = new PrintStream(filename.replace(':', '_').replace(' ', '_')+ "charge_" + c +  "_chargeOff_" + co + "_comp_" + com+".m");
							out.println("off=[");	
							for(InterPeakOffsetPeak p : offsetPeaksforOutput){
								if(p.gof.getBaseCharge() == c && p.gof.getChargeOffset() == co && p.gof.isComplementary() == com)
									out.println(p.gof.getOffset() + "\t" + p.y);
							}
							out.println("];");
							out.close();
						}
					}
				}
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		return offsetPeaks;
	}
	
	
	
	
		
	static public HashSet<InterPeakOffset> getGeneralizedOffsets(Peak bp, ArrayList<Peak> cps, boolean isComplementary, int chargeOffset, Tolerance tol){
		HashSet<InterPeakOffset> offs = new HashSet<InterPeakOffset>();
		for(Peak cp : cps){
			InterPeakOffset off = getGeneralizedOffset(bp, cp, isComplementary, chargeOffset, tol);
			if(off != null) offs.add(off);
		}
		
		return offs;
	}

	
	static public InterPeakOffset getGeneralizedOffset(Peak bp, Peak cp, boolean isComplementary, int chargeOffset, Tolerance tol){
		// base peaks already have charge offset	
		float offset = (cp.getMz() - bp.getMz()) * bp.getCharge();
	//	if(bp.getCharge() == 2) System.out.println(offset);
		if(offset > MAX || offset < MIN) return null;
		return new InterPeakOffset(offset, isComplementary, chargeOffset, bp.getCharge() - chargeOffset, tol.getToleranceAsDa(500)*2);
	}
	
}
