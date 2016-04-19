package tdaAnalysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class PostProcessOutput {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		boolean isISB = false;
		boolean xTandem = false;
		
		
		String[] proteinKeys1 = {"P", "Q", "A", "Y"};
		String[] proteinKeys2 = {"Y", "A", "I"};
		
		String[] proteinKeys = proteinKeys1;
		if(!isISB) proteinKeys = proteinKeys2;
		
		int pepCol = 7;
		int proteinCol = 8;
		int snCol = 1;
		
		if(xTandem){
			pepCol = 2;
			proteinCol = 3;
			snCol = 1;
		}
		
		String folder = "/home/kwj/workspace/inputs/TDA/NewKyowon" + (xTandem? "_X" : "") +   "/";
		if(isISB) folder = "/home/kwj/workspace/inputs/TDA/Kyowon" + (xTandem? "_X" : "") +   "/";
		
	
		for(int i=25;i<=27;i++){
			
			int maxSN = 4966;
			if(i == 22 || i == 23 || i == 28) maxSN = 47292;
			
			if(!isISB){
				maxSN = 9758;
				if(i == 22 || i == 23|| i == 28) maxSN = 55391;
			}
			
			//if(i == 27) maxSN = Integer.MAX_VALUE;
			
			if(!new File(folder + i).exists()) continue;
		
			for(File outFile : new File(folder + i).listFiles()){
				if(!xTandem && !outFile.getName().endsWith("txt"))continue;
				if(xTandem && !outFile.getName().endsWith("tsv"))continue;
				
			
				for(File dbFile : new File(folder + i + "/db").listFiles()){
					if(!dbFile.getName().endsWith("fasta"))continue;
					if(dbFile.getName().endsWith("reverse.fasta") || dbFile.getName().endsWith("Rsep.fasta")) continue;
					
					String db = dbFile.getAbsolutePath();
					String resultName = outFile.getAbsolutePath();
					if(resultName.endsWith("_p")) continue;
					
					BufferedLineReader in = new BufferedLineReader(resultName);
					PrintStream out = new PrintStream(resultName+"_p");
					SuffixArray sa = new SuffixArray(new SuffixArraySequence(db));
					
					String s;
					
					while((s=in.readLine())!=null){
						if(s.startsWith("#") || s.startsWith("Spec")){
							out.println(s);
							continue;
						}
						String[] token = s.split("\t");
							
					//	if(token[proteinCol].startsWith(proteinKeys[0])){
					//		out.println(s);
					//		continue;
					//	}
						String newProtein = null;
						
					//	if(i == 27 && Integer.parseInt(token[snCol]) > maxSN) continue;
						
						if(Integer.parseInt(token[snCol]) <= maxSN){

							String pep = token[pepCol];
							
							if(pep.contains("."))
								pep = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
						
						
								
							for(String proteinKey : proteinKeys){
								for(String protMatch : sa.getAllMatchingAnnotations(pep))
								{
									if(protMatch.startsWith(proteinKey)){
										newProtein = protMatch.split("\\s+")[0];
										break;
									}
								}
								if(newProtein != null) break;
							}
						}
						if(newProtein==null){
							out.println(s);
						}else{
							out.print(token[0]);
							for(int j=1; j<token.length; j++)
							{
								if(j != proteinCol)
									out.print("\t"+token[j]);
								else
									out.print("\t"+newProtein);
							}
							out.println();
							System.out.println("**************" + i + "\t" + s + "\t" + newProtein);
						}
						
						
					}
					
					
					in.close();
					out.close();
				}
				
				
			}
		}
		
		
		
		
		
	}

}
