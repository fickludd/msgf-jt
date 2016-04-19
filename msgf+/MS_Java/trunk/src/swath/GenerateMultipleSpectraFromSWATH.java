package swath;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import msutil.Peak;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MzXMLSpectraIterator;

public class GenerateMultipleSpectraFromSWATH {
	
	static private HashMap<Integer, ArrayList<Peak>> getPrecursors(String file) throws IOException{
		HashMap<Integer, ArrayList<Peak>> table = new HashMap<Integer, ArrayList<Peak>>();
		BufferedLineReader in = new BufferedLineReader(file);
		String s;
		int key = 0;
		ArrayList<Peak> value = null;
		
		while((s=in.readLine())!=null){
			String[] token = s.split("\t");
			if(s.startsWith("S")){
				if(value!=null) table.put(key, value);
				
				key = Integer.parseInt(token[1]);
				value = new ArrayList<Peak>();
			}else{
				int charge = Integer.parseInt(token[2]);
				float intensity = Float.parseFloat(token[3]);
				float mz = Float.parseFloat(token[4]);
				Peak p = new Peak(mz, intensity, charge);
				value.add(p);
			}	
		}
		table.put(key, value);
		
		
		in.close();
		
		return table;
	}
	
	static public void main(String[] args) throws IOException{
		String swath = "/home/kwj/workspace/inputs/SWATH/sw2/14304_EIF4A2_SWATH-1x.mzXML";
		String hardklor = "/home/kwj/workspace/inputs/SWATH/Ecoli/14304_EIF4A2_SWATH-1x_hardklor.hk";
		PrintStream out1 = new PrintStream("/home/kwj/workspace/inputs/SWATH/spectraFromSWATH/14304_EIF4A2_SWATH-1x_ver1.mgf");
		PrintStream out2 = new PrintStream("/home/kwj/workspace/inputs/SWATH/spectraFromSWATH/14304_EIF4A2_SWATH-1x_ver2.mgf");
		
		
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(swath);
		
		HashMap<Integer, ArrayList<Peak>> table = getPrecursors(hardklor);
		
		ArrayList<Integer> sns = new ArrayList<Integer>();
		
		sns.addAll(table.keySet());
		sns.add(Integer.MAX_VALUE);
		
		Collections.sort(sns);
		int i = 0;
		
		ArrayList<Peak> precursors = table.get(sns.get(i));
		
		while(itr.hasNext()){
			Spectrum s = itr.next();
			
			if(s.getScanNum() >= sns.get(i+1)){
				i++;
				precursors = table.get(sns.get(i));
			}
			
			if(precursors.isEmpty()) continue;
			
			ArrayList<Peak> sprecursors = new ArrayList<Peak>();
			
			for(Peak pp : precursors){
				float pm = (int)(s.getPrecursorPeak().getMz()/25)*25.0f;
				if(pp.getMz() >= pm && pp.getMz() < pm+25){
					sprecursors.add(pp);
					Peak newPrecursor = s.getPrecursorPeak().clone();
					newPrecursor.setCharge(pp.getCharge());
					newPrecursor.setMz(pp.getMz());
					
					s.setPrecursor(newPrecursor);
					
					s.outputMgf(out1);
				}
			}
			
			if(!sprecursors.isEmpty()){
				s.outputMgf(out2, sprecursors);
			}
			
		}
		
		out1.close();
		out2.close();
		
	}
	
}
