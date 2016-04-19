package kyowon_local;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Peptide;
import parser.BufferedLineReader;

public class gettol {
	public static void main(String[] args) throws IOException {
		
		BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/ROC/msdictionarySearchResults/Decoy/merged_new.txt");
		
		BufferedWriter out = new BufferedWriter(new FileWriter("/home/kwj/workspace/outputs/ROC/msdictionarySearchResults/Decoy/merged_new_0.5da.txt"));
		//BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/ROC/msdictionarySearchResults/merged_new.txt");
		//BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/ROC/inspectSearchResults/merged_pvalued.out");
		//BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/ROC/omssaSearchResults/whole.csv");
		
		
		String s;
		int pepcol = 2;
		int mzcol = 4;
		int chargecol = 5;
		float max = -100;
		String dil = "\t";
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#") || s.startsWith("Spectrum")) continue;
			String[] token = s.split(dil);
			
			if(Integer.parseInt(token[chargecol]) != 2) continue;
			
			float mz = Float.parseFloat(token[mzcol]);
			String p = (String) token[pepcol].subSequence(token[pepcol].indexOf(".")+1, token[pepcol].lastIndexOf("."));
			Peptide pep = new Peptide(p, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
			
			float mass = (mz - (float)Composition.PROTON) * (float)2; 
			//System.out.println(Math.abs(mass - pep.getParentMass())/2 + "\t" + token[mzcol+1]);
			//if(Math.abs(mass - pep.getParentMass()) < 0.5)
			//	out.write(s+"\n");
			
			max = Math.max(max, Math.abs(mass - pep.getParentMass()));
		}
		
		System.out.println(max);
		in.close();
		out.close();

	}

}
