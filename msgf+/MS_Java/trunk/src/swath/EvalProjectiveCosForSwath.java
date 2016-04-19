package swath;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import msutil.Peak;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MzXMLSpectraIterator;

public class EvalProjectiveCosForSwath {

	
	public static void main(String[] args) throws IOException{
		String swfile = "/home/kwj/workspace/inputs/SWATH/sw2/ACG_EIF4a2_swath_1x_vs_IDAlibrary_msplit.txt";
		String diafile = "/home/kwj/workspace/inputs/SWATH/sw2/ACG_msplitsearch_vs_IDAspeclib_projectedcosine.txt";
		String filekey = "sec_exclusion";//"no_exclusion";//"no_exclusion";//"";//"";// "";//"";//;
		String ioncntfile = "/home/kwj/workspace/inputs/SWATH/new/ACG_msplitsearch_vs_IDAspeclib.txt";
		
		String out = "/home/kwj/Dropbox/SWATH/report3/CosEval/" + filekey + "coseval.m";
		ArrayList<HashMap<Integer, ArrayList<Float>>> mSwath =  getlibQScoreMap(swfile, null,  null, null, true);
		ArrayList<HashMap<Integer, ArrayList<Float>>> mDIA =  getlibQScoreMap(diafile, null, "/home/kwj/workspace/inputs/SWATH/new/", filekey, true);
		PrintStream outp = new PrintStream(out);
		outp.println("a=[");
		for(int sn : mSwath.get(0).keySet()){
			if(mDIA.get(0).containsKey(sn)){
				outp.println(mSwath.get(0).get(sn).get(0) + "\t" + mDIA.get(0).get(sn).get(0) + "\t" + mSwath.get(1).get(sn).get(0) + "\t" + mDIA.get(1).get(sn).get(0) + "\t"+  mSwath.get(2).get(sn).get(0) + "\t" + mDIA.get(2).get(sn).get(0) + "\t" + sn);
			}
		}
		outp.println("];");
		outp.close();
	}
	
	public static ArrayList<HashMap<Integer, ArrayList<Float>>>  getlibQScoreMap(String file, String ioncntFile, String specFolder, String filekey, boolean selectBestLib) throws IOException {
		

		HashMap<String, ArrayList<Integer>> qLibMap = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, ArrayList<Float>> qLibScoreMap = new HashMap<String, ArrayList<Float>>();
		HashMap<String, ArrayList<Float>> qIonCntMap = new HashMap<String, ArrayList<Float>>();
		HashMap<String, ArrayList<Float>> qExpIonMap = new HashMap<String, ArrayList<Float>>();
		
		HashMap<Integer, ArrayList<String>> libQMap = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<Float>> libQScoreMap = new HashMap<Integer, ArrayList<Float>>(); // final product
		HashMap<Integer, ArrayList<Float>> libIonCntMap = new HashMap<Integer, ArrayList<Float>>();
		HashMap<Integer, ArrayList<Float>> libExpIonMap = new HashMap<Integer, ArrayList<Float>>();
		
		
		HashMap<String, ArrayList<Float>> ionCntMap = new HashMap<String, ArrayList<Float>>();
		
		if(specFolder!=null){
			for(File sf : new File(specFolder).listFiles()){
				if(!sf.getName().endsWith("mzXML")) continue;
				if(filekey!=null && !sf.getName().contains(filekey)) continue; 
				MzXMLSpectraIterator itr = new MzXMLSpectraIterator(sf.getAbsolutePath());
				
				while(itr.hasNext()){
					Spectrum sp = itr.next();					
					
					String sk = String.valueOf(sp.getPrecursorPeak().getMz());

					ionCntMap.put(sk, new ArrayList<Float>());
					
					float tp = 0;
					for(Peak pp : sp)
						tp += pp.getIntensity();
					ionCntMap.get(sk).add(tp);					
				}
			}
		}//14184_EIF4A2_A_fixedCE_TOF56k_P94_msgfdbIDs.mgf:836.36505
		String s;
		BufferedLineReader in;
		
		if(ioncntFile!=null){
			in = new BufferedLineReader(ioncntFile);
			while((s=in.readLine())!= null){
				if(s.startsWith("#") || s.startsWith("Parsed") || s.contains("DECOY")) continue;
				String[] token = s.split("\t");
				if(token.length < 5) continue;
				
				if(filekey!=null && !token[0].contains(filekey)) continue;
				
				float sk = Float.parseFloat(token[1]);
				if(ionCntMap.containsKey(sk))
					ionCntMap.get(sk).add(Float.parseFloat(token[token.length-1]));					
				
			}
			
			in.close();
		}
		in = new BufferedLineReader(file);
		
		
		while((s=in.readLine())!= null){
			if(s.startsWith("#") || s.startsWith("Parsed") || s.contains("DECOY")) continue;
				
			String[] token = s.split("\t");
			if(token.length < 5) continue;
			
			String[] t = token[0].split("/");
			
			String key = t[t.length-1] + ":" +	token[1];
			
			if(filekey!=null && !key.contains(filekey)) continue;
			
			int libSN = Integer.parseInt(token[8].split(":")[1].trim());
			float score = Float.parseFloat(token[7]); 
			float icn = Float.parseFloat(token[token.length-1]);
			//ionCntMap
			
			if(score > 0.99 || score < 0.6) continue; 
			
			float expIon = 0;
			if(specFolder != null){
				String sk = String.valueOf(Float.parseFloat(token[2]));
				//0.8247162	0.7579607	530645.0	7064.0	0.0	1726.0812	26274
				//icn = ionCntMap.get(sk).get(1);
				expIon = ionCntMap.get(sk).get(0) * Float.parseFloat(token[token.length-3]);
			}
			
			
			if(!qLibMap.containsKey(key)){
				qLibMap.put(key, new ArrayList<Integer>());
				qLibScoreMap.put(key, new ArrayList<Float>());
				qIonCntMap.put(key, new ArrayList<Float>());
				qExpIonMap.put(key, new ArrayList<Float>());
			}
			
			qLibMap.get(key).add(libSN);
			qLibScoreMap.get(key).add(score);
			qIonCntMap.get(key).add(icn);
			qExpIonMap.get(key).add(expIon);
		}
		
		if(selectBestLib){
			for(String key : qLibScoreMap.keySet()){
				float maxScore = 0;
				int index = 0;
				for(int i=0; i<qLibScoreMap.get(key).size();i++){
					float score= qLibScoreMap.get(key).get(i);
					if(score > maxScore){
						maxScore = score;
						index = i;
					}
				}
				
				ArrayList<Float> u = new ArrayList<Float>();
				u.add(maxScore);
				qLibScoreMap.put(key, u);
				
				ArrayList<Float> w = new ArrayList<Float>();
				w.add(qIonCntMap.get(key).get(index));
				qIonCntMap.put(key, w);
				
				ArrayList<Float> h = new ArrayList<Float>();
				h.add(qExpIonMap.get(key).get(index));
				qExpIonMap.put(key, h);
				
				ArrayList<Integer> v = new ArrayList<Integer>();
				v.add(qLibMap.get(key).get(index));
				qLibMap.put(key, v);
			}
		}
		
		for(String key : qLibMap.keySet()){
			for(int i=0;i<qLibMap.get(key).size();i++){
				int sn = qLibMap.get(key).get(i);
				if(!libQMap.containsKey(sn)){
					libQMap.put(sn, new ArrayList<String>());
					libQScoreMap.put(sn, new ArrayList<Float>());
					libIonCntMap.put(sn, new ArrayList<Float>());
					libExpIonMap.put(sn, new ArrayList<Float>());
				}
				libQMap.get(sn).add(key);
				libQScoreMap.get(sn).add(qLibScoreMap.get(key).get(i));
				libIonCntMap.get(sn).add(qIonCntMap.get(key).get(i));
				libExpIonMap.get(sn).add(qExpIonMap.get(key).get(i));
				
			}
		}
		
		for(int key : libQScoreMap.keySet()){
			float maxScore = 0;
			int index = 0;
			for(int i=0; i<libQScoreMap.get(key).size();i++){
				float score= libQScoreMap.get(key).get(i);
				if(score > maxScore){
					maxScore = score;
					index = i;
				}
			}
			
			ArrayList<Float> u = new ArrayList<Float>();
			u.add(maxScore);
			libQScoreMap.put(key, u);
			
			ArrayList<Float> w = new ArrayList<Float>();
			w.add(libIonCntMap.get(key).get(index));
			libIonCntMap.put(key, w);
			
			ArrayList<Float> h = new ArrayList<Float>();
			h.add(libExpIonMap.get(key).get(index));
			libExpIonMap.put(key, h);
			
			ArrayList<String> v = new ArrayList<String>();
			v.add(libQMap.get(key).get(index));
			libQMap.put(key, v);
		}
		
		in.close();
		
		ArrayList<HashMap<Integer, ArrayList<Float>>> ret = new ArrayList<HashMap<Integer, ArrayList<Float>>>();
		ret.add(libQScoreMap);
		ret.add(libIonCntMap);     
		ret.add(libExpIonMap);
		return ret;
	}	

}