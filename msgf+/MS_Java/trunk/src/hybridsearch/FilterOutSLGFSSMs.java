package hybridsearch;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;

import parser.BufferedLineReader;
import parser.PSM;
import parser.PSMList;
import parser.SLGFParser;

public class FilterOutSLGFSSMs {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String slgf = args[0];
		String msgfT = args[1];
		String out = args[2];
		
		HashSet<String> ssmKeys = new HashSet<String>();
		
		
		PSMList<PSM> slgfPSMs = SLGFParser.parse(slgf);
		for(PSM p : slgfPSMs){
			ssmKeys.add(p.getSpecFileName()+"\t"+p.getScanNum());
		}
		
		String s;
		PrintStream outs = new PrintStream(out);
		
		BufferedLineReader in = new BufferedLineReader(msgfT);
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#") || s.startsWith("Scan#")){
				outs.println(s);
				continue;
			}
			String[] token = s.split("\t");
			String k =  token[0].substring(1+token[0].lastIndexOf(System.getProperty("file.separator"))) + "\t" + token[1];
			if(ssmKeys.contains(k)) continue;
			outs.println(s);
			
		}
		
		in.close();
		outs.close();
		
	}

}
