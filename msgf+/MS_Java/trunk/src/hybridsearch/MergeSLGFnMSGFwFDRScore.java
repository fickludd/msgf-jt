package hybridsearch;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;

import parser.MSGFDBParser;
import parser.PSM;
import parser.PSMList;
import parser.SLGFParser;
import fdr.TargetDecoyPSMSet;

public class MergeSLGFnMSGFwFDRScore {

	private static TreeMap<Float,Float> getFDRMap(PSMList<PSM> list, boolean isGreaterBetter){

		ArrayList<Float> targetScores = new ArrayList<Float>();
		ArrayList<Float> decoyScores = new ArrayList<Float>();
		
		for(PSM p : list){
			if(!Float.isNaN(p.getScore("DECOY"))){
				if(p.getScore("DECOY")==1){
					decoyScores.add(p.getProbScore());
				}else{
					targetScores.add(p.getProbScore());
				}
			}else{
				if(p.getProtein().startsWith("REV")){
					decoyScores.add(p.getProbScore());
				}else{
					targetScores.add(p.getProbScore());
				}
			}
		}
		
		return TargetDecoyPSMSet.getFDRMap(targetScores, decoyScores, isGreaterBetter, true, 1);
		
	}
	

	public static void main(String[] args) throws IOException {
		String slgft = args[0];
		String slgfd = args[1];
		String msgft = args[2];
		String msgfd = args[3];
		
		PrintStream outt = new PrintStream(args[4]);
		PrintStream outd = new PrintStream(args[5]);
		
		PSMList<PSM> slgfTargetPSMs = SLGFParser.parse(slgft);
		PSMList<PSM> slgfDecoyPSMs = SLGFParser.parse(slgfd);
		PSMList<PSM> msgfTargetPSMs = MSGFDBParser.parse(msgft);
		PSMList<PSM> msgfDecoyPSMs = MSGFDBParser.parse(msgfd);
		
		PSMList<PSM> slgfPSMs = new PSMList<PSM>();
		slgfPSMs.addAll(slgfTargetPSMs);
		slgfPSMs.addAll(slgfDecoyPSMs);
		
		for(PSM p : slgfPSMs){
			p.probScore(-p.getProbScore());
		}
		
		PSMList<PSM> slgfSpecPSMs = slgfPSMs.getDistinctiveSpectralSet();
		PSMList<PSM> slgfPepPSMs = slgfSpecPSMs.getDistinctivePeptideSet();
		
		for(PSM p : slgfSpecPSMs){
			p.probScore(-p.getProbScore());
		}
		
		for(PSM p : slgfPepPSMs){
			p.probScore(-p.getProbScore());
		}
		
		
		TreeMap<Float,Float> slgfSpecFDRMap = getFDRMap(slgfSpecPSMs, true);
		TreeMap<Float,Float> slgfPepFDRMap = getFDRMap(slgfPepPSMs, true);
		
		
		PSMList<PSM> msgfPSMs = new PSMList<PSM>();
		msgfPSMs.addAll(msgfTargetPSMs);
		msgfPSMs.addAll(msgfDecoyPSMs);
		
		PSMList<PSM> msgfSpecPSMs = msgfPSMs.getDistinctiveSpectralSet();
		PSMList<PSM> nsgfPepPSMs = msgfSpecPSMs.getDistinctivePeptideSet();
		
		TreeMap<Float,Float> msgfSpecFDRMap = getFDRMap(msgfSpecPSMs, false);
		TreeMap<Float,Float> msgfPepFDRMap = getFDRMap(nsgfPepPSMs, false);
		
		
		for(PSM p : slgfTargetPSMs){
			float prob = p.getProbScore();
			float fdr = slgfSpecFDRMap.lowerEntry(prob).getValue();
			float pepFdr = slgfPepFDRMap.lowerEntry(prob).getValue();
			
			p.probScore(fdr);
			p.score("P", pepFdr);
			p.score("V", prob);
		}
		
		
		for(PSM p : slgfDecoyPSMs){
			float prob = p.getProbScore();
			float fdr = slgfSpecFDRMap.lowerEntry(prob).getValue();
			float pepFdr = slgfPepFDRMap.lowerEntry(prob).getValue();
			
			p.probScore(fdr);
			p.score("P", pepFdr);
			p.score("V", prob);
		}
		
		for(PSM p : msgfTargetPSMs){
			float prob = p.getProbScore();
			float fdr = msgfSpecFDRMap.higherEntry(prob).getValue();
			float pepFdr = msgfPepFDRMap.higherEntry(prob).getValue();
			
			p.probScore(fdr);
			p.score("P", pepFdr);
			p.score("V", prob);
		}
		
		
		for(PSM p : msgfDecoyPSMs){
			float prob = p.getProbScore();
			float fdr = msgfSpecFDRMap.higherEntry(prob).getValue();
			float pepFdr = msgfPepFDRMap.higherEntry(prob).getValue();
			
			p.probScore(fdr);
			p.score("P", pepFdr);
			p.score("S", prob);
		}
		
		
		PSMList<PSM> allTarget = new PSMList<PSM>();
		allTarget.addAll(msgfTargetPSMs);
		allTarget.addAll(slgfTargetPSMs);
		
		PSMList<PSM> allDecoy = new PSMList<PSM>();
		allDecoy.addAll(msgfDecoyPSMs);
		allDecoy.addAll(slgfDecoyPSMs);
		
		allTarget = allTarget.getDistinctiveSpectralSet();
		allDecoy = allDecoy.getDistinctiveSpectralSet();
		
		Collections.sort(allTarget, new PSM.PSMProbScoreComparator());
		Collections.sort(allDecoy, new PSM.PSMProbScoreComparator());
		
		outt.println("#SpecFile\tScan#\tPrecursor\tCharge\tPeptide\tisDecoy\tMSGFSpecProb\tSLGFScore\tHybridScore\tFDRScore\tPepFDRScore");
		
		for(PSM p : allTarget){
			if(Float.isNaN(p.getScore("V"))){ // msgf
				outt.println(p.getSpecFileName()+"\t"+p.getScanNum()+"\t"+p.getPrecursorMz()+"\t"+p.getCharge()+"\t"+
						p.getPeptideStr() + "\t" + p.getProtein().startsWith("REV") + "\t" + p.getScore("S") +"\tNaN\tNaN\t" + p.getProbScore() +  "\t" + p.getScore("P"));
			}else{
				outt.println(p.getSpecFileName()+"\t"+p.getScanNum()+"\t"+p.getPrecursorMz()+"\t"+p.getCharge()+"\t"+
						p.getPeptideStr() + "\t" + (p.getScore("DECOY") == 1) + "\tNaN\t" + p.getScore("V") +"\tNaN\t"  + p.getProbScore() +  "\t" + p.getScore("P"));
		
			}
		}

		outd.println("#SpecFile\tScan#\tPrecursor\tCharge\tPeptide\tisDecoy\tMSGFSpecProb\tSLGFScore\tHybridScore\tFDRScore\tPepFDRScore");
		
		for(PSM p : allDecoy){
			if(Float.isNaN(p.getScore("V"))){ // msgf
				outd.println(p.getSpecFileName()+"\t"+p.getScanNum()+"\t"+p.getPrecursorMz()+"\t"+p.getCharge()+"\t"+
						p.getPeptideStr() + "\t" + p.getProtein().startsWith("REV") + "\t" + p.getScore("S") +"\tNaN\tNaN\t" + p.getProbScore() +  "\t" + p.getScore("P"));
			}else{
				outd.println(p.getSpecFileName()+"\t"+p.getScanNum()+"\t"+p.getPrecursorMz()+"\t"+p.getCharge()+"\t"+
						p.getPeptideStr() + "\t" + (p.getScore("DECOY") == 1) + "\tNaN\t" + p.getScore("V") +"\tNaN\t"  + p.getProbScore() +  "\t" + p.getScore("P"));
		
			}
		}
		
		outt.close();
		outd.close();
		
		

	}

}
