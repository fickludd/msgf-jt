package swath;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;

public class SpecProbIonCntHisto {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String file = "/home/kwj/workspace/inputs/SWATH/new/ACG_msplitsearch_vs_IDAspeclib.txt";
		String file2 = "/home/kwj/workspace/inputs/SWATH/new/MSGFDB-757fe6ce-group_by_spectrum-main.txt";
		HashMap<String, Float> icMap = new HashMap<String, Float>();
		HashMap<String, Integer> pepNumTable = new HashMap<String, Integer>();
		HashMap<String, ArrayList<ArrayList<Float>>> ionCosTable = new HashMap<String, ArrayList<ArrayList<Float>>>();
		String out  = "/home/kwj/workspace/inputs/SWATH/new/ACG_msplitsearch_vs_IDAspeclib2.m";
		
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(file);
			String s;
			PrintStream ps = new PrintStream(out);
			ps.println("c=[");
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				
				String[] token = s.split("\t");
				
				//if(Float.parseFloat(token[8]) < 0.3f) continue;
				String[] t = token[0].split("/");
				
				
				String key = t[t.length-1] + ":" +	token[1];
				
				icMap.put(key, Float.parseFloat(token[token.length-1]));
				
			}
			
			in.close();
			
			in = new BufferedLineReader(file2);
	
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");
				String key = token[0] + ":" + token[1];
				if(!icMap.containsKey(key)) continue;
				
				float ic = icMap.get(key);
				
				String key2 = token[6] + ":" + token[7].substring(token[7].indexOf(".")+1, token[7].lastIndexOf("."));
				
				
				if(!ionCosTable.containsKey(key2)){
					ionCosTable.put(key2, new ArrayList<ArrayList<Float>>());
					pepNumTable.put(key2, 0);
				}
				
				ArrayList<Float> v = new ArrayList<Float>();
				//.out.println(s);
				

				v.add((float)Math.log10(ic));
				v.add((float)-Math.log10(Float.parseFloat(token[11])));
				ps.println(v.get(0) + "\t" + v.get(1));
				ionCosTable.get(key2).add(v);
				pepNumTable.put(key2, pepNumTable.get(key2)+1);
			}
			
			String k=null;
			int max = 0;
			
			for(String t : pepNumTable.keySet()){
				int n = pepNumTable.get(t);
				if(n > max && t.equals("3:FWDM+15.995AQLRAVVVDDYR")){//GIDVQQVSLVINYDLPTNR
					max = n;
					k = t;
				}
			}
			
			System.out.println("a=[");
			for(ArrayList<Float> nn : ionCosTable.get(k)){
				System.out.println(nn.get(0) + "\t" + nn.get(1));
			}
			System.out.println("];\n" + k);
			ps.println("];");
			ps.close();
			in.close();
			
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
	}

}
