package kyowon_local;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import parser.BufferedLineReader;

public class inspect {
	public static void main(String[] argv) throws IOException{
		BufferedLineReader in = new BufferedLineReader(argv[0]);
		BufferedWriter out = new BufferedWriter(new FileWriter(argv[0]+".target"));
		
		BufferedWriter out_d = new BufferedWriter(new FileWriter(argv[0]+".decoy"));
		
		String s;
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#")){
				out.write(s+"\n");
				out_d.write(s+"\n");
				continue;
			}
			
			String[] token = s.split("\t");
			if(Integer.parseInt(token[4]) != 2) continue;
			if(token[3].contains("XXX"))
				out_d.write(s+"\n");
			else
				out.write(s+"\n");
		}
		
		out.close();
		out_d.close();
		
		in.close();
		
		
	}
}
