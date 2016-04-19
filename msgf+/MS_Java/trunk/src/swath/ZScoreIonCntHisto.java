package swath;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msutil.Peak;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MzXMLSpectraIterator;

public class ZScoreIonCntHisto {

	static private float getMean(ArrayList<Float> xs){
		float m = 0;
		for(float x : xs) m+=x;
		
		return m/xs.size();
	}
	
	static private float getSD(ArrayList<Float> xs, float m){
		float v = 0;
		float ms = m*m;
		
		for(float x : xs) v+=x*x - ms;
		
		v /= (xs.size()-1);
		return (float)Math.sqrt(v);
	}
	
	
	static private float getAvgZscore(Spectrum s, float noiseFraction){
		ArrayList<Float> noise = new ArrayList<Float>();
		ArrayList<Float> signal = new ArrayList<Float>();
		int r = (int)(s.size()*(1-noiseFraction));
		
		s.setRanksOfPeaks();
		for(Peak p : s){
			if(p.getRank() < r) signal.add(p.getIntensity());
			else noise.add(p.getIntensity());
		}
		
		float m = getMean(noise);
		float sd = getSD(noise, m);
		float z = 0;
		
		for(float x : signal){
			z += (x-m)/sd;
		}
		
		return z/signal.size();
		
	}
	
	static private float getNumPeaksWithZscoreLargerThan(Spectrum s, float zScoreThreshold, float noiseFraction){
		ArrayList<Float> peaks = new ArrayList<Float>();
		ArrayList<Float> noise = new ArrayList<Float>();
		
		int r = (int)(s.size()*(1-noiseFraction));
		
		s.setRanksOfPeaks();
		for(Peak p : s){
			if(p.getRank() > r)  noise.add(p.getIntensity());
			peaks.add(p.getIntensity());
		}
		
		float m = getMean(noise);
		float sd = getSD(noise, m);
		float z = 0;
		
		for(float x : peaks){
			if((x-m)/sd>=zScoreThreshold) z++ ;
		}
		
		return z;
		
	}
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		float noiseFraction = .5f;
		float zscoreThreshold = 10;
		
		String xmlFolder = "/home/kwj/workspace/inputs/SWATH/sw2/";
		String file = xmlFolder + "ACG_EIF4a2_swath_1x_vs_IDAlibrary_msplit.txt";
	
		String file3 = "/home/kwj/workspace/inputs/SWATH/new/MSGFDB-757fe6ce-group_by_spectrum-main.txt";
		
		HashMap<String, Float> icMap = new HashMap<String, Float>();
		HashMap<String, ArrayList<Float>> expIonIntensityMap = new HashMap<String, ArrayList<Float>>();
		HashMap<String, ArrayList<Float>> cosMap = new HashMap<String, ArrayList<Float>>();
		
		HashSet<String> sns = new HashSet<String>();
	//	HashMap<String, Integer> pepNumTable = new HashMap<String, Integer>();
	//	HashMap<String, ArrayList<ArrayList<Float>>> ionZscoreTable = new HashMap<String, ArrayList<ArrayList<Float>>>();
		
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(file);
			String s;
			
	
			while((s=in.readLine())!=null){
				if(s.startsWith("#") || s.startsWith("Parsed")) continue;
				
				String[] token = s.split("\t");

				if(token.length < 5)continue;
				
				//if(Float.parseFloat(token[8]) < 0.3f) continue;
				String[] t = token[0].split("/");
				
				String key = t[t.length-1] + ":" +	token[1];
				icMap.put(key, Float.parseFloat(token[token.length-1]));
				
				
				if(!cosMap.containsKey(key)){
					cosMap.put(key, new ArrayList<Float>());
					expIonIntensityMap.put(key, new ArrayList<Float>());
				}
				expIonIntensityMap.get(key).add(Float.parseFloat(token[12]));
				cosMap.get(key).add(Float.parseFloat(token[7]));
				if(Float.parseFloat(token[7])<0.6f) System.out.println("Hmm");  
			}
			
			in.close();
			
			
			/*in = new BufferedLineReader(file3);
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");
				sns.add(token[0]+":"+token[1]);
			}
			*/
			int i = 0;
			for(File f : new File(xmlFolder).listFiles()){
				if(f.getName().endsWith("mzXML")){
					
					String out  = xmlFolder + "swathfiltered"+i+"z.m";
					PrintStream ps = new PrintStream(out);
					ps.println("'"+f.getName().replace("_", "")+"'");
					ps.println("z" + i + "=[");
					Iterator<Spectrum> iterator = new MzXMLSpectraIterator(f.getAbsolutePath());
					while(iterator.hasNext()){
						Spectrum spec= iterator.next();
						//
						
						float pi  = 0;
						for(Peak p: spec) pi += p.getIntensity();
						
						String key = f.getName() + ":" + spec.getScanNum();
					//	if(!sns.contains(key)) continue;
						
						if(icMap.containsKey(key)){
							//ps.println((float)Math.log10(icMap.get(key)) + "\t" + getNumPeaksWithZscoreLargerThan(spec, zscoreThreshold, noiseFraction) + "\t" + spec.getScanNum());
							for(int l=0; l<cosMap.get(key).size(); l++){
								float c  = cosMap.get(key).get(l);
								float expIon = expIonIntensityMap.get(key).get(l);
								ps.println((float)Math.log10(expIon * pi)  + "\t" + c + "\t" + spec.getScanNum() + "\t" + (float)Math.log10(icMap.get(key)));
							}
						}
						
					}
					ps.println("];");				
					ps.close();
					i++;
				}
			}
			
			
			in.close();
			
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
	}


}


