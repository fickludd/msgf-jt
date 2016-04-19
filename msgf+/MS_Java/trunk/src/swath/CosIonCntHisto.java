package swath;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;

public class CosIonCntHisto {
	static public void main(String[] args){
		HashMap<String, Integer> pepNumTable = new HashMap<String, Integer>();
		HashMap<String, ArrayList<ArrayList<Float>>> ionCosTable = new HashMap<String, ArrayList<ArrayList<Float>>>();
		
		String file = "/home/kwj/workspace/inputs/SWATH/sw2/ACG_swathdevelopment_IDA_vs_IDA_library_plusDecoy_projcos.txt";
		String out  = "/home/kwj/Dropbox/SWATH/hist2/ACG_swathdevelopment_IDA_vs_IDA_library_plusDecoy_projcos.m";
		String filter = "Control";
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(file);
			String s;
			PrintStream ps = new PrintStream(out);
			ps.println("b=[");
			while((s=in.readLine())!=null){
				if(s.startsWith("#") || s.startsWith("Parsed") || s.contains("DECOY")) continue;
				
				
				String[] token = s.split("\t");
				
				
				if(token.length < 5) continue;
				
				if(filter!=null && !token[0].contains(filter)) continue;
				
				if(Float.parseFloat(token[7]) > 0.99f) continue;
				
				String key = token[6]+":"+token[4];
				
				if(!ionCosTable.containsKey(key)){
					ionCosTable.put(key, new ArrayList<ArrayList<Float>>());
					pepNumTable.put(key, 0);
				}
				
				ArrayList<Float> v = new ArrayList<Float>();
				//.out.println(s);
				
				v.add((float)Math.log10(Float.parseFloat(token[token.length-1])));
				v.add(Float.parseFloat(token[7]));
				ps.println(v.get(0) + "\t" + v.get(1));
				ionCosTable.get(key).add(v);
				pepNumTable.put(key, pepNumTable.get(key)+1);
				
			}
			
			String k=null;
			int max = 0;
			
			for(String t : pepNumTable.keySet()){
				int n = pepNumTable.get(t);
				if(n > max && !t.equals("3:GIDVQQVSLVINYDLPTNR")){
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
			in.close();
			ps.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
}
