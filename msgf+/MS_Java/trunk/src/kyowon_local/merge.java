package kyowon_local;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import parser.BufferedLineReader;

public class merge {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File folder = new File(args[0]);
		File[] files = folder.listFiles();
		HashMap<Integer, ArrayList<String>> map = new HashMap<Integer, ArrayList<String>>();
		BufferedWriter out = new BufferedWriter(new FileWriter(args[1]));
		
		int maxsn = Integer.MIN_VALUE;
		///data/kyjeong/outputs/inspectSearchResults
		for(File file: files){
			if(!file.getName().endsWith(".txt")) continue;
			
			String s;
			BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
			while((s=in.readLine())!=null){
				if(s.startsWith("#"))continue;
				
				String[] token = s.split("\t");
				
				if(Float.parseFloat(token[14]) < 8) continue;
				
				int sn = Integer.parseInt(token[1]);
				
				ArrayList<String> v;
				if(map.containsKey(sn)){
					v = map.get(sn);
				}else
					v = new ArrayList<String>();
				
				v.add(s);
				
				map.put(sn, v);
				
				maxsn = Math.max(maxsn, sn);
			}
			in.close();
		}
		
		for(int sn=0; sn<=maxsn; sn++){
			if(map.containsKey(sn)){
				for(String s : map.get(sn)){
					out.write(s+"\n");
				}
			}
		}
		
		
		
		out.close();
		

	}

}
