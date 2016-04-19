package swath;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import msutil.Peak;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MzXMLSpectraIterator;

public class ExpIonCntvsTotalIonCnt {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		//pep -> abundance, best scan number ,rt, exp int, frac
		String swath = "/home/kwj/workspace/inputs/SWATH/Ecoli/ACG_swathdevelopment_14343_UPS1-400fm_Ecolilysate_vs_IDA_library_projcos_min0.6_pepFDR_filtered.txt";
		String mzXML = "/home/kwj/workspace/inputs/SWATH/Ecoli/spectra/" +
				"14343_UPS1_400fm_Ecolilysate_IDA_5600.mzXML";
		PrintStream out1 = new PrintStream(swath+"_table.txt");
		PrintStream out2 = new PrintStream(swath+"_tableMatlab.txt");
		
		HashMap<String, String> pepMap = new HashMap<String, String>();
		HashMap<Integer, ArrayList<Float>> specMap = new HashMap<Integer, ArrayList<Float>>();
		
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(mzXML);
		
		while(itr.hasNext()){
			Spectrum sp = itr.next();
			ArrayList<Float> v = new ArrayList<Float>();
			v.add(sp.getRt());
			
			float totalInt = 0;
			for(Peak p : sp) totalInt += p.getIntensity();
			
			v.add(totalInt);
			specMap.put(sp.getScanNum(), v);
			
		}
		System.out.println("Step 1");
		
		BufferedLineReader in = new BufferedLineReader(swath);
		
		String s;
		while((s=in.readLine())!=null){
			if(s.startsWith("#") || s.startsWith("Parsed") || s.startsWith("matching") || s.contains("DECOY")) continue;
			String[] token = s.split("\t");
			
			if(!token[token.length-1].contains("HUMAN")) continue;
			if(token[token.length-1].toLowerCase().contains("keratin")) continue;
			
			
			String pep = token[4];
			float cos = Float.parseFloat(token[7]);
			//if(cos < 0.81f) continue;
			//if(Float.parseFloat(token[11]) <= 15) continue;
			
			if(pepMap.containsKey(pep)){
				String[] btoken = pepMap.get(pep).split("\t");
				if(cos < Float.parseFloat(btoken[0])) continue;
			}
			
			String[] subToken = token[token.length-1].split(" ");
			float abundance = Float.parseFloat(subToken[subToken.length-1]);
			
			int sn = Integer.parseInt(token[1]);
			float frac = Float.parseFloat(token[12]);
			
			ArrayList<Float> v = specMap.get(sn);
			float rt = v.get(0);
			float totalInt = v.get(1);
			
			float totalIcn = Float.parseFloat(token[token.length-2]);
			
			String value = cos + "\t" + abundance + "\t" + sn + "\t" + rt + "\t" + totalIcn + "\t" +  totalInt + "\t" + (totalInt*frac) + "\t" + frac;
				//pep -> abundance, best scan number ,rt, exp int, frac
			pepMap.put(pep, value);
		}
		
		out1.println("#Peptide\tProjected Cosine Score\tAbundance\tScan#\tRetention Time\tIon Count\tTotal Peak Intensity\tExplained Peak Intensity\tFraction of Explained Peak Intensity");
		out2.println("a=[");
		
		for(String pep : pepMap.keySet()){
			String v = pepMap.get(pep);
			String[] token = v.split("\t");
			out1.println(pep+"\t"+v);
			out2.println(token[6] + "\t" + token[7]);
			
		}
		out2.println("];");
		
		in.close();
		out1.close();
		out2.close();
	}

}
