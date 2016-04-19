package hybridsearch;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeMap;

import parser.PSM;
import parser.PSMList;
import parser.SLGFParser;
import fdr.TargetDecoyPSMSet;

public class LibraryDBSearch {

	/**
	 * @param args
	 */
	private static void printUsage(){
		System.out.println("java -jar -Xmx3g LibraryDBSearch.jar [SLGFTarget] [SLGFDecoy] [MSGFTarget] [MSGFDecoy] [outputFileNamePrefix] [FDR] [FDR option] [Combine option]");
		
		System.out.println("SLGFTarget : output file from SLGF search against target only");
		System.out.println("SLGFDecoy : output file from SLGF search against decoy only");
		System.out.println("MSGFTarget : output file from MSGF search against target only");
		System.out.println("MSGFDecoy : output file from MSGF search against decoy only");
		
		System.out.println("outputFileNamePrefix : if [Combine option] is 0, two files [outputFileNamePrefix+ \"SLGF.txt\"] and [outputFileNamePrefix+\"MSGF.txt\"] will be generated. Each contains IDs from each tool.");
		System.out.println("                       else if [Combine option] is 1, one file [outputFileNamePrefix+ \".txt\"] will be generated.");
		
		System.out.println("FDR : 0-1 (only for combine option 0)");
		System.out.println("FDR option : 0 - multiple PSMs (SSMs) per spectrum, 1 - best PSM (SSM) per spectrum, 2 - peptide level FDR");
		System.out.println("Combine option : 0 - combine method 0, 1 - combine method 1");
		System.out.println();
		System.exit(0);
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length != 8) printUsage();
		
		String SLGFT = args[0];
		String SLGFD = args[1];
		String MSGFT = args[2];
		String MSGFD = args[3];
		
		String outputPrefix = args[4];
		
		float fdr = Float.parseFloat(args[5]);
		int fdrOption = Integer.parseInt(args[6]);
		int combineOption = Integer.parseInt(args[7]);
			
		PSMList<PSM> targetSLGF = SLGFParser.parse(SLGFT);
		PSMList<PSM> decoySLGF = SLGFParser.parse(SLGFD);
	
		PSMList<PSM> targetMSGF = SLGFParser.parse(MSGFT);
		PSMList<PSM> decoyMSGF = SLGFParser.parse(MSGFD);
		
		for(PSM psm : targetSLGF){
			psm.probScore(-psm.getProbScore()); // stupid but should keep..
			psm.score("T/D", 0);
		}
		for(PSM psm : decoySLGF){
			psm.probScore(-psm.getProbScore());
			psm.score("T/D", 1);
		}
		for(PSM psm : targetMSGF){			
			psm.score("T/D", 0);
		}
		for(PSM psm : decoyMSGF){			
			psm.score("T/D", 1);
		}
		
		PSMList<PSM> aSLGF = new PSMList<PSM>();
		
		aSLGF.addAll(targetSLGF);
		aSLGF.addAll(decoySLGF);
					
		PSMList<PSM> aMSGF = new PSMList<PSM>();
		
		aMSGF.addAll(targetMSGF);
		aMSGF.addAll(decoyMSGF);
		
		if(fdrOption >= 1){
			aSLGF = aSLGF.getDistinctiveSpectralSet();
			aMSGF = aMSGF.getDistinctiveSpectralSet();
		}
		if(fdrOption == 2){
			aSLGF = aSLGF.getDistinctivePeptideSet();
			aMSGF = aMSGF.getDistinctivePeptideSet();
		}
		
		
		if(combineOption == 0){
			PrintStream outSLGF = new PrintStream(outputPrefix+"SLGF.txt");
			PrintStream outMSGF = new PrintStream(outputPrefix+"MSGF.txt");
			
			HashSet<Integer> identifiedBySLGF = new HashSet<Integer>();
			
			PSMList<PSM> tSLGF = new PSMList<PSM>();
			PSMList<PSM> dSLGF = new PSMList<PSM>();
			
			for(PSM psm : aSLGF){
				if(psm.getScore("T/D") == 1) dSLGF.add(psm);
				else tSLGF.add(psm);
			}
			
			PSMList<PSM> tMSGF = new PSMList<PSM>();
			PSMList<PSM> dMSGF = new PSMList<PSM>();
						
			for(PSM psm : aMSGF){
				if(identifiedBySLGF.contains(psm.getScanNum())) continue;
				if(psm.getScore("T/D") == 1) dMSGF.add(psm);
				else tMSGF.add(psm);
			}
			
			
			PSMList<PSM> idSLGF = PSMList.selectUsingFDR(tSLGF, dSLGF, fdr);
			
			System.out.println("# identified SLGF PSMs at FDR "+fdr*100+"% : "+ idSLGF.size());
			
			for(PSM p : idSLGF){
				identifiedBySLGF.add(p.getScanNum());
			}
			
			
			
			PSMList<PSM> idMSGF = PSMList.selectUsingFDR(tMSGF, dMSGF, fdr);
			System.out.println("# identified MSGF PSMs at FDR "+fdr*100+"% : "+ idMSGF.size());
			
			outSLGF.println("#SpecFile\tScan#\tPrecursor\tCharge\tPeptide\tpValueScore");
			for(PSM p : idSLGF){
				outSLGF.println(p.getSpecFileName()+"\t"+p.getScanNum()+"\t"+p.getPrecursorMz()+"\t"+p.getCharge()
						+"\t"+p.getPeptideStr()+"\t"+(-p.getProbScore()));
			}
			
			
			outMSGF.println("#SpecFile\tScan#\tPrecursor\tCharge\tPeptide\tProtein\tSpecProb");
			for(PSM p : idMSGF){
				outMSGF.println(p.getSpecFileName()+"\t"+p.getScanNum() + "\t" + p.getPrecursorMz() + "\t" + p.getCharge()
						+ "\t" + p.getPeptideStr() + "\t" + p.getProtein() + "\t" + p.getProbScore());
			}
			
			outSLGF.close();
			outMSGF.close();
			
		}else {
			PrintStream out = new PrintStream(outputPrefix+".txt");

			ArrayList<Float> tscoreSLGF = new ArrayList<Float>(); 
			ArrayList<Float> dscoreSLGF = new ArrayList<Float>();
			ArrayList<Float> tscoreMSGF = new ArrayList<Float>();
			ArrayList<Float> dscoreMSGF = new ArrayList<Float>();
			
			for(PSM psm : aSLGF){
				if(psm.getScore("T/D") == 1) dscoreSLGF.add(psm.getProbScore());
				else tscoreSLGF.add(psm.getProbScore());
			}
								
			for(PSM psm : aMSGF){
				if(psm.getScore("T/D") == 1) dscoreMSGF.add(psm.getProbScore());
				else tscoreMSGF.add(psm.getProbScore());
			}
			
			TreeMap<Float,Float> slgfFDRMap = TargetDecoyPSMSet.getFDRMap(tscoreSLGF, dscoreSLGF, false, true, 1);
			TreeMap<Float,Float> msgfFDRMap = TargetDecoyPSMSet.getFDRMap(tscoreMSGF, dscoreMSGF, false, true, 1);
			
			
			for(PSM psm : targetSLGF){
				float fdrScore = slgfFDRMap.higherEntry(psm.getProbScore()).getValue();
				psm.probScore(fdrScore); 
			}
			for(PSM psm : decoySLGF){
				float fdrScore = slgfFDRMap.higherEntry(psm.getProbScore()).getValue();
				psm.probScore(fdrScore); 
			}
			for(PSM psm : targetMSGF){			
				float fdrScore = msgfFDRMap.higherEntry(psm.getProbScore()).getValue();
				psm.probScore(fdrScore);
			}
			for(PSM psm : decoyMSGF){			
				float fdrScore = msgfFDRMap.higherEntry(psm.getProbScore()).getValue();
				psm.probScore(fdrScore); 
			}
			
			
			PSMList<PSM> all = new PSMList<PSM>();
			all.addAll(targetSLGF);all.addAll(decoySLGF);all.addAll(targetMSGF);all.addAll(decoyMSGF);
			
			if(fdrOption >= 1){
				all = all.getDistinctiveSpectralSet();
			}
			if(fdrOption == 2){
				all = all.getDistinctivePeptideSet();
			}
			
			out.println("#SpecFile\tScan#\tPrecursor\tCharge\tPeptide\tempirical FDR");
			for(PSM p : all){
				out.println(p.getSpecFileName()+"\t"+p.getScanNum() + "\t" + p.getPrecursorMz() + "\t" + p.getCharge()
						+ "\t" + p.getPeptideStr() + "\t" + p.getProbScore());
			}
			
			out.close();
		}

	}

}
