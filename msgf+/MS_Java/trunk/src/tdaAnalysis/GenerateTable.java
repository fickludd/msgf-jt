package tdaAnalysis;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

import msgf.Tolerance;
import parser.MSGFDBParser;
import parser.PSM;
import parser.PSMList;
import parser.XTandemParser;

public class GenerateTable {
	
	static private boolean[] getFilter(int searchNum, boolean isISB){
		boolean[] ret = new boolean[3];
		
		if(searchNum != 1 && searchNum != 2 && searchNum != 9 && searchNum != 14 && searchNum != 16
				 && searchNum != 19  && searchNum != 23){
			ret[0] = true;
		}
		
		if(isISB && searchNum != 16 && searchNum != 17 && searchNum != 18){
		//if(searchNum == 1 || searchNum == 2||searchNum == 9 ||searchNum == 16 ||searchNum == 19 ){// searchNum != 17 && searchNum != 18){
			ret[1] = true;
		}
		
		
		if(searchNum != 1 && searchNum != 2 && searchNum != 3 && searchNum != 4 && searchNum != 9 && searchNum != 12
				 && searchNum != 16 && searchNum != 19){
			ret[2] = true;
		}
		
		return ret;
	}
	
	
	private FactualFDR getScoreThreshold(PSMList<PSM> thits, PSMList<PSM> dhits, PSMList<PSM> falseHitsAmongPassedHits,  String folder, String effFolder, int searchNum, boolean isISB, boolean isMSGFDB, boolean pepLevelFDR, float ffdrThreshold, int charge, boolean fixFFDR, String keyword){
		
		boolean[] filters = getFilter(searchNum, isISB);
		//String specKey = "mix"; 
		
		int maxSN = 4966;
		if(searchNum == 22 || searchNum == 23) maxSN = 47292;
		
		if(!isISB){
		//	specKey = "CPTAC";
			maxSN = 9758;
			if(searchNum == 22 || searchNum == 23) maxSN = 55391;
		}
		
		
		String[] proteinKeys;
		if(isISB){
			proteinKeys = new String[2];
			proteinKeys[0] = "P";
			proteinKeys[1] = "Q";
		//	proteinKeys[2] = "REV_P";
		//	proteinKeys[3] = "REV_Q";
		}else{
			proteinKeys = new String[1];
			proteinKeys[0] = "Y";
		//	proteinKeys[1] = "REV_Y";
		}
		
		boolean protFilter = filters[0]; 
		boolean pmFilter = filters[1];
		boolean specFilter = filters[2];
		
		
		
		PSMList<PSM> result = new PSMList<PSM>();
		PSMList<PSM> effResult = new PSMList<PSM>();
		
		FactualFDR ff = null;
		
		
		for(File file : new File(folder).listFiles()){
			if(!file.getName().endsWith("_p")) continue;
			
			if(!keyword.isEmpty() && !file.getName().contains(keyword)) continue;
			
			if(isMSGFDB) result.addAll(MSGFDBParser.parse(file.getAbsolutePath()));
			else result.addAll(XTandemParser.parse(file.getAbsolutePath()));
			
		}
		
		if(searchNum !=9 && searchNum !=10 && searchNum !=11) result = result.getDistinctiveSpectralSet();
		
		//System.out.println(result.size());
		
		
		for(File file : new File(effFolder).listFiles()){
			if(searchNum == 22 || searchNum == 23){
				if(!file.getName().endsWith("_q"))continue;
			}else{
				if(!file.getName().endsWith("_p")) continue;
			}
			
			if(!keyword.isEmpty() && !file.getName().contains(keyword)) continue;
			if(keyword.equals("filtered") && file.getName().contains("unfiltered")) continue; 
			if(keyword.isEmpty() && file.getName().contains("filtered")) continue;
			
			//if(file.getName().contains("filtered")) continue;
			
			
			if(isMSGFDB) 
				effResult.addAll(MSGFDBParser.parse(file.getAbsolutePath()));
			else effResult.addAll(XTandemParser.parse(file.getAbsolutePath()));
			
		}
		
		
		effResult = effResult.getDistinctiveSpectralSet();
		
		float ffdr = 1, pepffdr = 1, fdr = 1, pepfdr = 1;
		float t = 1e-5f;
		
		if(!isMSGFDB) t = -0f;
		
		if(isMSGFDB){
			if(searchNum == 19 ||searchNum == 20 ||searchNum == 21) t = -10f;
		}else if(searchNum == 19 ||searchNum == 20 ||searchNum == 21) t = 10f;
		
		PSMList<PSM> ehits = new PSMList<PSM>();
		
		//FactualFDR 
		ff = null;
		
		PSMList<PSM> tmpThits = new PSMList<PSM>(thits);
		PSMList<PSM> tmpDhits = new PSMList<PSM>(dhits);
		PSMList<PSM> tmpFalseHitsAmongPassedHits = new PSMList<PSM>(falseHitsAmongPassedHits);
		
		while(pepLevelFDR ? (fixFFDR? pepffdr > ffdrThreshold : pepfdr > ffdrThreshold) : (fixFFDR? ffdr > ffdrThreshold: fdr > ffdrThreshold)){
			ehits.clear();
			
			tmpThits = new PSMList<PSM>(thits);
			tmpDhits = new PSMList<PSM>(dhits);
			tmpFalseHitsAmongPassedHits = new PSMList<PSM>(falseHitsAmongPassedHits);
			//tmpThits.clear(); tmpDhits.clear();  tmpFalseHitsAmongPassedHits.clear();
			
			if(isMSGFDB){
				if(searchNum == 19 ||searchNum == 20 ||searchNum == 21){
					t += 0.1;
				}else{
					t/=1.05f;
				}
			}else{
				if(searchNum == 19 ||searchNum == 20 ||searchNum == 21){
					t/=1.05f;
				}else{
					t += 0.1;
				}
			}
			
			
			
			
			//HashSet<Integer> sns = new HashSet<Integer>();
			
			for(PSM psm : result){
				
				if(charge > 0 && charge < 4 && psm.getCharge() != charge) continue;
				if(charge >= 4 &&  psm.getCharge() < 4) continue;
				
				//sns.add(psm.getScanNum());
				
				if(isMSGFDB){
					if(searchNum == 19 ||searchNum == 20 ||searchNum == 21){
						if(psm.getRawScore() < t) continue;
					}else{
						if(psm.getProbScore() > t) continue;
					}
				}else{
					if(searchNum == 19 ||searchNum == 20 ||searchNum == 21){
						if(psm.getProbScore() > t) continue;
					}else{
						if(psm.getRawScore() < t) continue;
					}
				}
				
				if(psm.getProtein().startsWith("R") || psm.getProtein().startsWith("S")){
					tmpDhits.add(psm);					
				}else{ 
					tmpThits.add(psm);					
				}
				
				
			}
			
			for(PSM psm : effResult){
				if(charge > 0 && charge < 4 && psm.getCharge() != charge) continue;
				if(charge >= 4 &&  psm.getCharge() < 4) continue;
				//if(!sns.contains(psm.getScanNum())) continue;
				
				
				if(isMSGFDB){
					if(searchNum == 19 ||searchNum == 20 ||searchNum == 21){
						if(psm.getRawScore() < t) continue;
					}else{
						if(psm.getProbScore() > t) continue;
					}
				}else{
					if(searchNum == 19 ||searchNum == 20 ||searchNum == 21){
						if(psm.getProbScore() > t) continue;
					}else{
						if(psm.getRawScore() < t) continue;
					}
				}
				
				if(psm.getProtein().startsWith("R") || psm.getProtein().startsWith("S")){
					ehits.add(psm);	
				}
				
				
				
				
				
			}
			
			
			ff = new FactualFDR(tmpThits, isISB, isMSGFDB);
			
			if(protFilter){
				ff.setProteinFilter(proteinKeys);
			}
			if(pmFilter){
				ff.setPMFilter(new Tolerance(50f, true));//TODO
			}
			if(specFilter){
				//if(isMSGFDB){
				ff.setSpecFilter(maxSN);
				//}
				//else{
					//ff.setSpecFilter(specKey);
				//}
			}
			ff.calculate();
			
			//FactualFDR ff2 = new FactualFDR(ehits, isISB, isMSGFDB);
			
			
			
			int nthits = tmpThits.getDistinctiveSpectralSet().size();
			int ndhits = tmpDhits.getDistinctiveSpectralSet().size();
			
			fdr = (float)ndhits/nthits;
			if(searchNum == 12 || searchNum == 13) fdr = (float)2 * ndhits/(nthits + ndhits);
			
			pepfdr = (float) tmpDhits.getDistinctivePeptideSet().size()/tmpThits.getDistinctivePeptideSet().size();
				
			//lseHitsAmongPassedHits = new PSMList<PSM>();//P3
					
			//HashMap<Integer, PSM> ffm = new HashMap<Integer, PSM>();
			HashMap<Integer, PSM> ehitsm = new HashMap<Integer, PSM>();
			
			//for(PSM h : ff.getPassedPSMs()){
			//	ffm.put(h.getScanNum(), h);
			//}
			for(PSM h : ehits){
				//if(ffm.containsKey(h.getScanNum()))
				ehitsm.put(h.getScanNum(), h);
			}
			
			//System.out.println(ehitsm.size() + "\t" + ffm.size() + "\t" + ehits.size());
			
			for(int sn : ehitsm.keySet()){
				//if(ffm.containsKey(sn)){
					//PSM h = ffm.get(sn);
					PSM h2 = ehitsm.get(sn);//TODO
					//if(h == null || (h.getScanNum() == h2.getScanNum() && h.getRawScore() < h2.getRawScore()))
					//if(h!=null)
					tmpFalseHitsAmongPassedHits.add(h2);	
				//	}
				//}
			}
			
			tmpFalseHitsAmongPassedHits.addAll(ff.getFalsePSMs());
			//System.out.println(ndhits + "\t" + ff.getFalsePSMs().getDistinctiveSpectralSet().size() + "\t" + falseHitsAmongPassedHits.getDistinctiveSpectralSet().size());
			//System.out.println(nthits + "\t" + ff.getFalsePSMs().getDistinctiveSpectralSet().size()
			//		+"\t" + ff.getPassedPSMs().getDistinctiveSpectralSet().size());
			
			ffdr = tmpFalseHitsAmongPassedHits.getDistinctiveSpectralSet().size();
			ffdr /= nthits;
			
			pepffdr = tmpFalseHitsAmongPassedHits.getDistinctivePeptideSet().size();
			pepffdr /= tmpThits.getDistinctivePeptideSet().size();
				
		}
		
		System.out.println(tmpThits.getDistinctiveSpectralSet().size()-thits.getDistinctiveSpectralSet().size());
		System.out.println(tmpFalseHitsAmongPassedHits.getDistinctiveSpectralSet().size()-falseHitsAmongPassedHits.getDistinctiveSpectralSet().size());
		
		HashSet<Integer> sns = new HashSet<Integer>();
		for(PSM p : thits) sns.add(p.getScanNum());
			
		for(PSM p : tmpThits)
			if(!sns.contains(p.getScanNum())) thits.add(p);
		sns.clear();

		for(PSM p : dhits) sns.add(p.getScanNum());
		
		for(PSM p : tmpDhits)
			if(!sns.contains(p.getScanNum()))  dhits.add(p);
		sns.clear();
		
		for(PSM p : falseHitsAmongPassedHits) sns.add(p.getScanNum());
		
		for(PSM p : tmpFalseHitsAmongPassedHits)
			if(!sns.contains(p.getScanNum())) falseHitsAmongPassedHits.add(p);
		
		
		System.out.println("Score Threhsold: " + t);
		return ff;
		
	}
	

	private void outputResult(int searchNum, boolean isISB, boolean isMSGFDB, boolean pepLevelFDR, float ffdrThreshold, int charge, boolean fixFFDR){
		
		String folder = "/home/kwj/workspace/inputs/TDA/Kyowon" + (!isMSGFDB? "_X" : "") +   "/" + searchNum + "/";
		
		if(!isISB)
			folder = "/home/kwj/workspace/inputs/TDA/NewKyowon" + (!isMSGFDB? "_X" : "") +   "/" + searchNum + "/";
		
		String effFolder = "/home/kwj/workspace/inputs/TDA/Kyowon" + (!isMSGFDB? "_X" : "") +   "/" + 28 + "/";
		
		if(!isISB)
			effFolder = "/home/kwj/workspace/inputs/TDA/NewKyowon" + (!isMSGFDB? "_X" : "") +   "/" + 28 + "/";
		
	 
		 PSMList<PSM> thitsFinal = new PSMList<PSM>();
		 PSMList<PSM> dhitsFinal = new PSMList<PSM>();
		 PSMList<PSM> falseHitsAmongPassedHitsFinal = new PSMList<PSM>();
		
		 String[] keywords = new String[1];
		 keywords[0] = "";
		 
		 if(searchNum == 25 || searchNum ==  27){
			 keywords = new String[2];
			 keywords[0] = "before";
			 keywords[1] = "filtered";
		 }
		 
		 if(searchNum == 24 || searchNum ==  26){
			 keywords = new String[1];
			 keywords[0] = "unfiltered";
		 }
		 
		 for(int i=0; i<keywords.length; i++){
			 String keyword = keywords[i];
			// PSMList<PSM> thits = new PSMList<PSM>();
			// PSMList<PSM> dhits = new PSMList<PSM>();
			// PSMList<PSM> falseHitsAmongPassedHits = new PSMList<PSM>();
			 
			
			 
			 getScoreThreshold(thitsFinal, dhitsFinal, falseHitsAmongPassedHitsFinal, folder, effFolder, searchNum, isISB, isMSGFDB, pepLevelFDR, ffdrThreshold, charge,  (i==keywords.length-1? fixFFDR : false), keyword);
			
			 //System.out.println(thitsFinal.size());
			 
			 //thitsFinal.addAll(thits);
			// dhitsFinal.addAll(dhits);
			// falseHitsAmongPassedHitsFinal.addAll(falseHitsAmongPassedHits);
			
		 }
		 

		 
		 int nthits = thitsFinal.getDistinctiveSpectralSet().size();
		 int ndhits = dhitsFinal.getDistinctiveSpectralSet().size();
		 
		 float ffdr = falseHitsAmongPassedHitsFinal.getDistinctiveSpectralSet().size();
		 ffdr /= nthits;
		 
		 
		 float pepffdr = falseHitsAmongPassedHitsFinal.getDistinctivePeptideSet().size();
		 pepffdr /= thitsFinal.getDistinctivePeptideSet().size();

	//	 for(PSM psm : falseHitsAmongPassedHitsFinal.getDistinctivePeptideSet()){
			// System.out.println(p);
	//		 float pm = (psm.getPrecursorMz() - (float)Composition.PROTON) * psm.getCharge(); 
				
			// System.out.println(Math.abs(pm - psm.getPeptide().getParentMass()) <= new Tolerance(30f, true).getToleranceAsDa(pm));				
	//	 }
		 
		 float fdr = (float)ndhits/nthits;
		 if(searchNum == 12 || searchNum == 13) fdr = (float)2 * ndhits/(nthits + ndhits);
			
		 float pepfdr = (float) dhitsFinal.getDistinctivePeptideSet().size()/thitsFinal.getDistinctivePeptideSet().size();
		
		// System.out.println( dhitsFinal.getDistinctivePeptideSet());
		// System.out.println( falseHitsAmongPassedHitsFinal.getDistinctivePeptideSet());
		 
		 int p = nthits;
		 int q = ndhits;
		 int r = (falseHitsAmongPassedHitsFinal.getDistinctiveSpectralSet().size());
			
		 int p1 = thitsFinal.getDistinctivePeptideSet().size();
		 int q1 = dhitsFinal.getDistinctivePeptideSet().size();
		 int r1 =(falseHitsAmongPassedHitsFinal.getDistinctivePeptideSet().size());

		 FisherExact fe = new FisherExact(100000000);

		 System.out.println("FactFDR = " + String.format("%.3f", ffdr));
		 System.out.println("FDR = " + String.format("%.3f\t%d\t%d\t%d", fdr, nthits, ndhits , falseHitsAmongPassedHitsFinal.getDistinctiveSpectralSet().size()));
		 System.out.println("FactPepFDR = " + String.format("%.3f", pepffdr));
		 System.out.println("PepFDR = " + String.format("%.3f", pepfdr));
		
		 System.out.println("ID = " + thitsFinal.getDistinctiveSpectralSet().size() + "\t#PEP = "+ thitsFinal.getDistinctivePeptideSet().size() + "\t" + dhitsFinal.getDistinctivePeptideSet().size());
		
		 if(searchNum == 12 || searchNum==13){
			 System.out.println("pValue (PSM level) = " + String.format("%.2f%%", 100 * fe.getCumlativeP(p+q,2*q,p,r)));
			 System.out.println("pValue (Peptide level) = " + String.format("%.2f%%", 100 * fe.getCumlativeP(p1,q1,p1,r1)));
 
		 }else{
			 System.out.println("pValue (PSM level) = " + String.format("%.2f%%", 100 * fe.getCumlativeP(p,q,p,r)));
			 System.out.println("pValue (Peptide level) = " + String.format("%.2f%%", 100 * fe.getCumlativeP(p1,q1,p1,r1)));
			 System.out.println("pValue (Peptide-PSM level) = " + String.format("%.2f%%", 100 * fe.getCumlativeP(p,q,p1,r1)));
			// System.out.println("pValue (Peptide-PSM level) = " + String.format("%.2f%%", 100 * fe.getCumlativeP(3529,45,2574,25)));
						
		 }
		}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		boolean isISB = false;
		boolean fixFFDR = false;
		// exclude 1 2 9 16 19 for yeast...
		
		int searchNum = 26;//
		
		int charge = 0;//
		
		boolean pepLevelFDR = false;//
		
		float ffdrThreshold = 0.010f;
		if(isISB) ffdrThreshold = 0.05f;
		GenerateTable gt = new GenerateTable();
		System.out.println("Search Number: " + searchNum);
		for(int i=0;i<2;i++){
			System.out.println("_____________________");
			gt.outputResult(searchNum, isISB, i==0, pepLevelFDR, ffdrThreshold, charge, fixFFDR);
		}
		
	}

}
