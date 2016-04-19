package swath;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import parser.BufferedLineReader;

public class CountPepts {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO AuHashMap<K, V>rated method stub
		HashMap<String, HashSet<String>> peps = new HashMap<String, HashSet<String>>();
		
		BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/inputs/SWATH/Ecoli/MSGFDB-c5e6c1f6-group_by_spectrum-main.tsv");
		String s;
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#")) continue;
			
			String[] token = s.split("\t");
			if(Float.parseFloat(token[14]) > 0.01 ) continue;
			if(!peps.containsKey(token[0])) peps.put(token[0], new HashSet<String>());
			
			peps.get(token[0]).add(token[7].substring(token[7].indexOf('.')+1, token[7].lastIndexOf('.')));
			
		}
		
		in.close();
		
		for(String k : peps.keySet()){
			System.out.println(k + "\t" + peps.get(k).size());
		}
		
 	}

}
