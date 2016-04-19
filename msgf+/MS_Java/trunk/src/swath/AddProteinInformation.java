package swath;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import parser.BufferedLineReader;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class AddProteinInformation {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String fasta = "/home/kwj/workspace/inputs/SWATH/Ecoli/Ecoli_genome_plusUPS_withProtQuant.fasta";
		String dir = "/home/kwj/workspace/inputs/SWATH/Ecoli/";
		
		for(File f : new File(dir).listFiles()){
			String result = f.getAbsolutePath();
			System.out.println(result);
			if(!result.endsWith(".txt") || result.endsWith("proteinAdded.txt")) continue;
			PrintStream newResult = new PrintStream(result+".proteinAdded.txt");
			
			SuffixArray sa;
			sa = new SuffixArray(new SuffixArraySequence(fasta));
			
			BufferedLineReader in = new BufferedLineReader(result);
			String s;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("Parsed") || s.startsWith("matching")){
					newResult.println(s);
					continue;
				}
				
				if(s.contains("DECOY")){
					newResult.println(s+"\tDECOY");
					continue;
				}
				
				if(s.startsWith("#")){
					newResult.println(s + "\tProtein");
					continue;
				}
				String[] token = s.split("\t");
				
				if(token.length<5) System.out.println(s);
				ArrayList<String> matchedProtSeq = sa.getAllMatchingAnnotations(token[4]);
				
				if(matchedProtSeq.isEmpty()){
					newResult.println(s+"\tUnmatched");
					System.out.println("OH NO " + s);
					continue;
				}
				newResult.println(s+"\t"+matchedProtSeq.get(0));
			}
			
			in.close();
			newResult.close();
		}
		
		
	}

}
