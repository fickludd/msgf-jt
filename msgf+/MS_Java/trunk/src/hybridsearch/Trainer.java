package hybridsearch;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import msutil.Peptide;
import parser.MSGFDBParser;
import parser.MSGFLibParser;
import parser.PSM;
import parser.PSMList;
import parser.SLGFParser;

public class Trainer {
	private String poorQualLibFile, msgfDir, paraFile, poorQualLibmsgfFile, specDir, goodQualmsgfFile, fasta;
	
	public Trainer(String poorQualLibFile, String goodQualmsgfFile, String fasta, String msgfDir, String specDir, String paraFile){
		this.poorQualLibFile = poorQualLibFile;
		this.poorQualLibmsgfFile = poorQualLibFile + "_msgf";
		this.fasta = fasta;
		this.goodQualmsgfFile = goodQualmsgfFile;
		this.specDir = specDir;
		this.msgfDir = msgfDir;
		this.paraFile = paraFile;
	}
	
	public void run(){
		//runMSGF();
		trainConditionalProbabilities();
		trainAlpha();
	}
	/*6Training for conditional probabilities ... done
Training for alpha = 1.0
6973	2341
 ... done : # ID = 3350
1.1798744E-8
8144	1567
 ... done : # ID = 4204
-0.191493
8003	1730
 ... done : # ID = 4853	 # ep = 239
3.932091E-7
MSGF: 3350
SLGF: 4204
Hybrid: 4853
MSGF/SLGF: 2704
Hybrid/SLGF: 4062
MSGF/Hybrid: 3188
Common: 2694
51358	92852
3228	1931
3.1837342E-9
SLGF + MSGF ID : 4712
Peps : 3977	4125
3795
*/
	private void trainAlpha(){
	//	System.out.println(160.03065-AminoAcidSet.getStandardAminoAcidSet().getAminoAcid('C').getMass());
	//	System.out.println(new Peptide("QFGALM+15.994915NESQASCDK").getMass());
	//	System.out.println(new Peptide("C-17.027EDVFEYKDDSAIK").getMass());
		
		
		int maxNumID = 0;
		float alpha = 0;
		for(float a = 1f;a<1f+0.01;a+=0.1){
			System.out.println("Training for alpha = " + a);
		//	MSGFLib ml = new MSGFLib(paraFile + a , goodQualmsgfFile, fasta, poorQualLibFile,  paraFile, specDir);
		//	MSGFLib.setAlpha(a);
		//	ml.run();
			
			PSMList<PSM> hybridPSMs = MSGFLibParser.parse(paraFile + a);
			
			PSMList<PSM> target = new PSMList<PSM>();
			PSMList<PSM> decoy = new PSMList<PSM>();
			
			PSMList<PSM> originalPSMs = MSGFDBParser.parse(goodQualmsgfFile);
			
			for(PSM psm : originalPSMs.getDistinctiveSpectralSet()){
				if(psm.getProtein().startsWith("REV")) decoy.add(psm);
				else target.add(psm);
			}
			System.out.println(target.size() + "\t" + decoy.size());
			//target = target.getDistinctivePeptideSet();
			target = PSMList.selectUsingFDR(target, decoy, 0.01f);
			int numID = target.size();
			System.out.println(" ... done : # ID = " + numID);
			System.out.println(target.get(target.size()-1).getProbScore());
			
			//target.clear(); decoy.clear();
			
			PSMList<PSM> target2 = new PSMList<PSM>();
			PSMList<PSM> decoy2 = new PSMList<PSM>();
			
			PSMList<PSM> originalPSMs2 = SLGFParser.parse(poorQualLibFile);
			
			for(PSM psm : originalPSMs2){
				psm.probScore(-psm.getProbScore());
			}
			
			
			
			for(PSM psm : originalPSMs2.getDistinctiveSpectralSet()){
				//psm.probScore(-psm.getProbScore());
				if(psm.getScore("T/D") == 1  ) decoy2.add(psm);
				else target2.add(psm);
			}
			System.out.println(target2.size() + "\t" + decoy2.size());
			//target2 = target2.getDistinctivePeptideSet();
			target2 = PSMList.selectUsingFDR(target2, decoy2, 0.01f);
			numID = target2.size();
			System.out.println(" ... done : # ID = " + numID);
			System.out.println(target2.get(target2.size()-1).getProbScore());
			
			//target.clear(); decoy.clear();
			PSMList<PSM> target3 = new PSMList<PSM>();
			PSMList<PSM> decoy3 = new PSMList<PSM>();
			
			int ep = 0;
			for(PSM psm : hybridPSMs.getDistinctiveSpectralSet()){
				if(psm.getProtein().startsWith("REV")) decoy3.add(psm);
				else{
					target3.add(psm);				
				}
			}
			System.out.println(target3.size() + "\t" + decoy3.size());

			//target3 = target3.getDistinctivePeptideSet();
			
			target3 = PSMList.selectUsingFDR(target3, decoy3, 0.01f); // TODO;
			numID = target3.size();
			for(PSM psm : target3){
				if(psm.getScore("IsEPvalue")==1) ep++;
			}
			
			System.out.println(" ... done : # ID = " + numID + "\t # ep = " + ep);
			System.out.println(target3.get(target3.size()-1).getProbScore());
			
			HashSet<Integer> sns = new HashSet<Integer>();
			for(PSM psm : target){
				sns.add(psm.getScanNum());
			}
			HashSet<Integer> sns2 = new HashSet<Integer>();
			for(PSM psm : target2){
				sns2.add(psm.getScanNum());
			}
			HashSet<Integer> sns3 = new HashSet<Integer>();
			for(PSM psm : target3){
				sns3.add(psm.getScanNum());
			}
			
			int c = 0;
			
			HashSet<Integer> tsns = new HashSet<Integer>();
			for(int sn : sns){
				if(sns3.contains(sn)) c++; 
				
			}

			int c2 = 0;
			for(int sn : sns2){
				if(sns3.contains(sn)) c2++; 
				else tsns.add(sn);
			}
			
			int c3 = 0;
			int c4 = 0;
			for(int sn : sns){
				if(sns2.contains(sn)){
					c3++;
					if(sns3.contains(sn)) c4 ++;
				}
			}
			
			for(PSM psm : target2){
				if(tsns.contains(psm.getScanNum())){
				//	System.out.println(psm.getScanNum() + "\t" + psm.getPeptideStr() + "\t" + psm.getProbScore());
				}
			}
			
			
			System.out.println("MSGF: " + (sns.size()));
			System.out.println("SLGF: " + (sns2.size()));
			System.out.println("Hybrid: " + (sns3.size()));
			
			System.out.println("MSGF/SLGF: " + (c3 ));
			System.out.println("Hybrid/SLGF: " + (c2));
			System.out.println("MSGF/Hybrid: " + (c));
			
			
			System.out.println("Common: " + c4);
			
			
			target.clear(); decoy.clear();
			PSMList<PSM> filteredPSMs = new PSMList<PSM>();
			for(PSM psm : originalPSMs){
				if(sns2.contains(psm.getScanNum())) continue;
				filteredPSMs.add(psm);
			}
			System.out.println(filteredPSMs.size() + "\t" + originalPSMs.size());
			for(PSM psm : filteredPSMs.getDistinctiveSpectralSet()){
				if(psm.getProtein().startsWith("REV")) decoy.add(psm);
				else target.add(psm);
			}
			System.out.println(target.size() + "\t" + decoy.size());
			
			target = PSMList.selectUsingFDR(target, decoy, 0.01f);
			System.out.println(target.get(target.size()-1).getProbScore());
			
			numID = target.size();
			System.out.println("SLGF + MSGF ID : " + (sns2.size()  +  numID));

			HashSet<Peptide> peps = new HashSet<Peptide>();
			for(PSM psm : target.getDistinctiveSpectralSet().getDistinctivePeptideSet()){
				peps.add(psm.getPeptide());
			}
			for(PSM psm : target2.getDistinctiveSpectralSet().getDistinctivePeptideSet()){
				peps.add(psm.getPeptide());
			}
			System.out.println("Peps : " +peps.size() + "\t" + target3.getDistinctivePeptideSet().size());
			c=0;
			for(Peptide pep : peps){
				boolean com = false;
				for(PSM psm1 : target3.getDistinctiveSpectralSet().getDistinctivePeptideSet()){
					
					//if(psm1.getScanNum() == 928) System.out.println(psm1.getPeptideStr());
					if(pep.equals(psm1.getPeptide())){
						com = true;
						c ++;
					}
				}
			//	if(!com) System.out.println(pep.toStringWithModification(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys()));
			}
			
			System.out.println(c);
		//	//System.out.println((sns.size() -c ) + " " + c + " " +  (sns3.size() - c));
			//System.out.println((sns2.size() -c2 ) + " " + c2+ " " +  (sns3.size() - c2));
			//System.out.println((sns.size() -c3 ) + " " + c3+ " " +  (sns2.size() - c3));
			
			if(numID > maxNumID){
				alpha = a;
				maxNumID = numID;
			}
		}
		
		try {
			PrintWriter ps = new PrintWriter(new FileWriter(paraFile, true));
			ps.println("#Alpha\t"+alpha);
			ps.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void runMSGF(){
		String[] args = {
			"java","-jar", "-Xmx2000M", msgfDir + System.getProperty("file.separator") + "MSGF.jar" ,
			"-i", poorQualLibFile, "-d", specDir , "-o", poorQualLibmsgfFile, // ETD HCD?
		}; 
        Process process;
		try {
			process = new ProcessBuilder(args).start();
		
	        InputStream is = process.getInputStream();
	        InputStreamReader isr = new InputStreamReader(is);
	        BufferedReader br = new BufferedReader(isr);
	        String line;
	
	        System.out.printf("Output of running %s is:", 
	           Arrays.toString(args));
	
	        while ((line = br.readLine()) != null) {
	          System.out.println(line);
        }
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	private void trainConditionalProbabilities(){
		System.out.print("Training for conditional probabilities");
		PSMList<PSM> psms = SLGFParser.parse(poorQualLibmsgfFile);
		HashMap<Integer, ArrayList<ConditionalProbability>> histogram = new HashMap<Integer, ArrayList<ConditionalProbability>>();
		try {
			PrintWriter ps = new PrintWriter(new FileWriter(paraFile));
			
			for(PSM psm : psms){
				float specProb = psm.getRawScore();
				float slgfPvalue = psm.getProbScore();
				HybridPvalue hp = new HybridPvalue(specProb, slgfPvalue);
				
				//if(hp.getIntegerLibScore() < 0)
				//	System.out.println(slgfPvalue);
				ConditionalProbability cp = new ConditionalProbability(hp.getIntegerDBScore(), hp.getIntegerLibScore());
				
				if(!histogram.containsKey(hp.getIntegerDBScore()))
					histogram.put(hp.getIntegerDBScore(), new ArrayList<ConditionalProbability>());
				
				histogram.get(hp.getIntegerDBScore()).add(cp);
				
			}
			
			for(int i : histogram.keySet()){
				float[] libScoreDist = new float[HybridPvalue.maxScore];
				int t =  histogram.get(i).size();
				
				for(ConditionalProbability cp : histogram.get(i)){
					libScoreDist[cp.getLibScore()] ++;
				}
				
				for(int j=0;j<libScoreDist.length;j++){
					float v = libScoreDist[j] / t;
					ps.println(i+"\t"+j+"\t"+v);
				}
				
			}
		
			ps.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		System.out.println(" ... done");
		
	}
	
	static public void main(String[] args){
		String poorQualLibFile = "/home/kwj/workspace/inputs/slgf/20080315_CPTAC6_16_6D014.mzXML.search_results_decoy.txt";
		String goodQualmsgfFile = "/home/kwj/workspace/inputs/slgf/msgfdb.txt";
		String fasta = "/home/kwj/workspace/inputs/DataBases/yeast.revConcat.fasta";
//		String poorQualLibFile = "/Users/kwj/Documents/workspace/inputs/annotatedHeckWholeSum_CID1.txt";
	//	String goodQualmsgfFile = "/Users/kwj/Documents/workspace/inputs/annotatedHeckWholeSum_ETD1.txt";
		
		String msgfDir = "/home/kwj/workspace/MSGF";
		String specDir = "/home/kwj/workspace/inputs/slgf/";
		String paraFile = "/home/kwj/workspace/inputs/slgf/outPara.txt";
		
		Trainer trainer = new Trainer(poorQualLibFile, goodQualmsgfFile, fasta, msgfDir, specDir, paraFile);
		trainer.run();
	}
	
}
