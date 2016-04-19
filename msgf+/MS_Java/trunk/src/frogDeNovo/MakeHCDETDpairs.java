package frogDeNovo;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class MakeHCDETDpairs {
	static public void main(String[] args) throws IOException{
		String HCD = "/home/kwj/workspace/inputs/Nebedev/HCD_MS2_20110103_03_Frog_init2_1-1000.mgf";
		String ETD = "/home/kwj/workspace/inputs/Nebedev/ETD_MS2_20110103_03_Frog_init2_1-1000.mgf";
		
		PrintStream H = new PrintStream(HCD+".mgf");
		PrintStream E = new PrintStream(ETD+".mgf");
		
		HashMap<Float, ArrayList<Spectrum>> tb = new HashMap<Float, ArrayList<Spectrum>>();
		
		Iterator<Spectrum> iH = new SpectraIterator(HCD, new MgfSpectrumParser());
		Iterator<Spectrum> iE = new SpectraIterator(ETD, new MgfSpectrumParser());
		
		while(iH.hasNext()){
			Spectrum s = iH.next();
			float pm = s.getPrecursorPeak().getMz()+
			Float.parseFloat(s.getTitle().split(" ")[12]);
			
			if(!tb.containsKey(pm)){
				tb.put(pm, new ArrayList<Spectrum>());
				tb.get(pm).add(s);
			}
			
		}
		
		while(iE.hasNext()){
			Spectrum s = iE.next();
			float pm = s.getPrecursorPeak().getMz() +
			Float.parseFloat(s.getTitle().split(" ")[12]);
			if(!tb.containsKey(pm)){
				tb.put(pm, new ArrayList<Spectrum>());
			}
			if(tb.get(pm).size() == 1)
				tb.get(pm).add(s);
		}
		
		for(float pm : tb.keySet()){
			ArrayList<Spectrum> ss = tb.get(pm);
			if(ss.size() != 2){
				System.out.println(ss);
				continue;
			}
			ss.get(0).outputMgf(H);			
			ss.get(1).outputMgf(E);
		}
		
		H.close();
		E.close();
		
	}
}
