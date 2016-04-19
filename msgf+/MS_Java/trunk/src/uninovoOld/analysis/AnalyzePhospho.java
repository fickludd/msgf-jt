package uninovoOld.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;

public class AnalyzePhospho {
	public static void main(String[] args) throws IOException{
		float threshold = 1.0e-15f;
		
		HashSet<String> targetHits = new HashSet<String>();
		HashSet<String> targetPhosphoHits = new HashSet<String>();
		HashSet<String> decoyHits = new HashSet<String>();
		
		String[] keyWords = { //"+80",
				"Y+80", "S+80", "T+80",
				};
		
		String dir = "/home/kwj/Desktop/CIDETD/";
		String out = dir + "Merged";
		PrintStream p = new PrintStream(out);
		
		for(File file : new File(dir).listFiles()){
			if(file.getName().endsWith("mod")){
				BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
				String s;
				while((s=in.readLine())!= null){
					p.println(s);
				}
				in.close();
			}
		}
		
		p.close();
		
		BufferedLineReader in = new BufferedLineReader(out);
		String s;
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#")) continue;
			String[] token = s.split("\t");
			if(Float.parseFloat(token[6]) > threshold) continue;
			if(token[8].startsWith("REV_"))
				decoyHits.add(token[5]);
			else{
				targetHits.add(token[5]);
				for(String key : keyWords){
					if(token[5].contains(key)){
						targetPhosphoHits.add(token[5]);
						break;
					}
				}
			}
		}
		in.close();
		System.out.println(targetHits.size() + "\t" + decoyHits.size() + "\t" + targetPhosphoHits.size());
	}
}
