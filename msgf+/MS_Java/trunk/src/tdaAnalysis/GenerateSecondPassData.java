package tdaAnalysis;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import msdbsearch.ReverseDB;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MSGFDBParser;
import parser.MgfSpectrumParser;
import parser.PSM;
import parser.PSMList;
import parser.XTandemParser;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class GenerateSecondPassData {
	static public void main(String[] args) throws Exception{
		
		boolean xTandem = true;
		
		String output = "/home/kwj/workspace/inputs/TDA/NewKyowon" + (xTandem? "_X" : "") +   "/7/XTANDEM_MIN-cfdd3da2-group_by_spectrum-main.tsv_p";//search output
		String db = "/home/kwj/workspace/inputs/TDA/NewKyowon" + (xTandem? "_X" : "") +   "/24/db/YAR.fasta"; // put concatenated for mode 1. put sep for mode 2, 3
		String mgf = "/home/kwj/workspace/inputs/TDA/NewKyowon" + (xTandem? "_X" : "") +   "/24/mgf/YA.mgf";
	
		float scoreThreshold =  33.00011f;//7.3805834E-10f;//31.0f;//7.521394E-10f
		
		PSMList<PSM> list = null;
		
		if(!xTandem) list = MSGFDBParser.parse(output);
		else list = XTandemParser.parse(output);
		
		HashSet<String> toput = new HashSet<String>();
    	HashSet<String> toput2 = new HashSet<String>();
    	HashSet<Integer> sns = new HashSet<Integer>();
    	
		for(int mode = 1;mode <=3; mode ++){
		
			toput.clear(); toput2.clear();
			
			
			PrintStream out = new PrintStream(db+"_filtered_" + mode + ".fasta");
			PrintStream out2 = null;
			
			int numT=0, numD=0;
			
			if(mode ==3)
				out2 = new PrintStream(db+"_filtered_rev_" + mode + ".fasta");
			
			SuffixArray sa = new SuffixArray(new SuffixArraySequence(db));
			
			
	    
			for(PSM psm : list){
				//System.out.println(psm.getRawScore());
				if(xTandem && (psm.getRawScore() < scoreThreshold)) continue;
				if(!xTandem && (psm.getProbScore() > scoreThreshold)) continue;
				
				sns.add(psm.getScanNum()-1);
				String pepStr = psm.getPeptideStr();
				
				ArrayList<String> matchedEntries = sa.getAllMatchingAnnotations(pepStr);
		    	ArrayList<String> matchedProtSeq = sa.getAllMatchingEntries(pepStr);
		    	
		    	for(int i=0; i<matchedEntries.size(); i++){
		    		if(mode != 1 && matchedEntries.get(i).startsWith("REV")){
		    			if(mode == 3){
		    				toput2.add(">"+matchedEntries.get(i)+"\n"+matchedProtSeq.get(i)+"\n");
		    			}
		    			continue;
		    		}
		    		toput.add(">"+matchedEntries.get(i)+"\n"+matchedProtSeq.get(i)+"\n");
		    	}
		    	
		    	
			}
			for(String t : toput){
	    		out.print(t);
	    		if(t.contains(">REV")) numD++;
	    		else numT++;	
			}
			if(mode == 1) System.out.println(numT + "\t" + numD);
			out.close();
			if(mode == 3){
				for(String t : toput2)
		    		out2.print(t);
		    	
				
			}
				
			if(mode == 2){// generate rev and conc
				ReverseDB.reverseDB(db+"_filtered_" + mode + ".fasta", db+"_filtered_rev_" + mode + ".fasta");
			}
			
			if(mode == 3){
				toput.clear();
				sa = new SuffixArray(new SuffixArraySequence(db+"_filtered_rev_" + 2 + ".fasta"));
				ArrayList<String> matchedEntries = sa.getAllMatchingAnnotations("M");
		    	ArrayList<String> matchedProtSeq = sa.getAllMatchingEntries("M");
		    	
		    	for(int i=0; i<matchedEntries.size(); i++){
		    		toput.add(">"+matchedEntries.get(i)+"\n"+matchedProtSeq.get(i)+"\n");
		    	}
		    	
				int tofill = toput.size() - toput2.size();
			
				for(String t : toput){
					if(tofill <=0 ) break;
					if(toput2.contains(t)) continue;
					
					out2.print(t);
					tofill--;
				}
				
				
				out2.close();
			}
		}
		Iterator<Spectrum> iterator = new SpectraIterator(mgf, new MgfSpectrumParser());
		
		PrintStream mgfout = new PrintStream(mgf+"_filtered.mgf");
		
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			if(sns.contains(s.getScanNum())) continue;
			s.outputMgf(mgfout);
		}
		mgfout.close();
		
		
	}
}
