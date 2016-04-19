package kyowon_local;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;

public class order {
	static public void main(String[] argv) throws IOException{
		BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/ROC/msdictionarySearchResults/Decoy/merged.txt");
		BufferedWriter out = new BufferedWriter(new FileWriter("/home/kwj/workspace/outputs/ROC/msdictionarySearchResults/Decoy/merged_new.txt"));
	
		HashMap<Integer, ArrayList<String>> map = new HashMap<Integer, ArrayList<String>>();
		
		String s;
		int max = Integer.MIN_VALUE;
		
		while((s=in.readLine())!=null){
			String[] token = s.split("\t");
			int sn = Integer.parseInt(token[1]);
			
			ArrayList<String> v;
			
			if(map.containsKey(sn))
				v = map.get(sn);
			else 
				v= new ArrayList<String>();
			
			v.add(s);
			
			map.put(sn, v);
			max = Math.max(max, sn);
		}
		in.close();
		
		for(int sn = 0; sn<=max;sn++){
			if(map.containsKey(sn)){
				for(String t : map.get(sn)){
					out.write(t+"\n");
				}
			}
		}
		out.close();
	}
}
