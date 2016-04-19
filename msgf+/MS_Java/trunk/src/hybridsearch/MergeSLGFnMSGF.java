package hybridsearch;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Collections;

import parser.MSGFDBParser;
import parser.PSM;
import parser.PSMList;
import parser.SLGFParser;

public class MergeSLGFnMSGF {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String slgf = args[0];
		String msgf = args[1];
		String out = args[2];
		
		PrintStream outs = new PrintStream(out);
		
		outs.println("#SpecFile\tScan#\tPrecursor\tCharge\tPeptide\tisDecoy\tMSGFSpecProb\tSLGFScore\tHybridScore\tFDR\tPepFDR");
		
		PSMList<PSM> slgfPSMs = SLGFParser.parse(slgf, true);
		PSMList<PSM> msgfPSMs = MSGFDBParser.parse(msgf, true);
		
		for(PSM psm : slgfPSMs){
			float f = psm.getScore("F");
			float p = psm.getProbScore();
			psm.probScore(f);
			psm.score("V", p); //  for sorting
		}
		
		for(PSM psm : msgfPSMs){
			float f = psm.getScore("F");
			float p = psm.getProbScore();
			psm.probScore(f);
			psm.score("S", p); //  for sorting
		}
		

		PSMList<PSM> total = new PSMList<PSM>();
		
		total.addAll(msgfPSMs);
		total.addAll(slgfPSMs);
		
		Collections.sort(total, new PSM.PSMProbScoreComparator());
		for(PSM p : total){
			if(Float.isNaN(p.getScore("V"))){ // msgf
				outs.println(p.getSpecFileName()+"\t"+p.getScanNum()+"\t"+p.getPrecursorMz()+"\t"+p.getCharge()+"\t"+
						p.getPeptideStr() + "\t" + p.getProtein().startsWith("REV") + "\t" + p.getScore("S") +"\tNaN\tNaN\t" + p.getProbScore() +  "\t" + p.getScore("P"));
			}else{
				outs.println(p.getSpecFileName()+"\t"+p.getScanNum()+"\t"+p.getPrecursorMz()+"\t"+p.getCharge()+"\t"+
						p.getPeptideStr() + "\t" + (p.getScore("DECOY") == 1) + "\tNaN\t" + p.getScore("V") +"\tNaN\t"  + p.getProbScore() +  "\t" + p.getScore("P"));
		
			}
		}

	
		
		
		
		outs.close();
		
		
	}

}
