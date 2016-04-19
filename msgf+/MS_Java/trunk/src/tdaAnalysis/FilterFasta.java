package tdaAnalysis;

import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class FilterFasta {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/inputs/TDA/NewKyowon_X/24/db/YAR.fasta_filtered_3.fasta");
		PrintStream out = new PrintStream("/home/kwj/workspace/inputs/TDA/NewKyowon_X/28/db/YR.fasta_filtered_3.fasta");
		
		String s;
		boolean towrite = true;
		while((s=in.readLine())!=null){
			if(s.startsWith(">A") || s.startsWith(">REV_A") || s.startsWith(">Y")){
				towrite = false;
			}else if(s.startsWith(">REV_Y")) towrite = true;
			
			if(towrite){
				out.println(s);
			}
			
		}
		
		out.close();
		in.close();
		
	}	
	

}
