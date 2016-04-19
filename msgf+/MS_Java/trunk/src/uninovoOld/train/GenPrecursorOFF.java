package uninovoOld.train;

import java.io.IOException;
import java.util.Iterator;

import msutil.Composition;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;

public class GenPrecursorOFF {
	
	static int range = 120;
	static int charge = 3;
	static WindowFilter filter = new WindowFilter(6, 50);
	
	static float getChargeChangedPrecursorMz(Peak precursor, int charge){
		return (precursor.getMass())/charge + (float)Composition.PROTON;
	}
	
	
	static public void main(String[] args) throws IOException{
		String specfilename = "/home/kwj/workspace/inputs/Training/ETD_Tryp_Confident.mgf";
		
		Iterator<Spectrum> iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
		
		float[][] offsets = new float[charge+1][range * 2+1];
		int specNum = 0;
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(spec.getCharge() != charge) continue;
			
			spec = filter.apply(spec);
			
			for(int c = 1; c <= charge ; c++){
				float mz = getChargeChangedPrecursorMz(spec.getPrecursorPeak(), c);
				for(Peak p : spec.getPeakListByMassRange(mz - range/c, mz + range/c)){ // mass offset
					float offset = (p.getMz() - mz) *c;
					offsets[c][Math.round(offset + range)]++;
				}
			}
			specNum++;
		}
		
		for(int c = 1; c <= charge; c++){
			System.out.println("off"+c+"=[");
			for(int i = 0; i<offsets[c].length; i++)
				System.out.println(i-range + "\t" + offsets[c][i]/specNum);
			System.out.println("];");
		}
		
	}
	
}
