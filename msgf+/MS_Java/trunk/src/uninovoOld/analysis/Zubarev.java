package uninovoOld.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;
import sequences.Constants;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;
import uninovoOld.DenovoReconstruction;
import uninovoOld.IPDGenerator;
import uninovoOld.SpectrumGraph;
import uninovoOld.UniNovo;

public class Zubarev {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		
		String specDir = "/home/kwj/workspace/inputs/Zubarev/";
		String fasta = "/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.v3.52.target.fasta";
		
		File[] specFiles = new File(specDir).listFiles();
		int mode = Integer.parseInt(args[0]);
		int maxGapNum = Integer.parseInt(args[1]);
		int numPerLength = Integer.parseInt(args[2]);		
		int minLength = Integer.parseInt(args[3]);
		
		int maxPTMNum = 1;
		String outputFile = "/home/kwj/Dropbox/HCDETDnew_"+mode + "_" + maxGapNum + "_" + numPerLength +  "_" + minLength + ".txt";
		if(args[4].startsWith("decoy")){
			fasta = "/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.v3.52.decoy.fasta";
			outputFile = "/home/kwj/Dropbox/HCDETDnew_"+mode + "_" + maxGapNum + "_" + numPerLength +  "_" + minLength + "_decoy.txt";
			
		}
		
		PrintStream op = new PrintStream(outputFile);
		
		float setAccuracyThreshold = 0.8f;
		int maxLength = 30;
		
		int minNumPerLength = 0;
		
		ArrayList<String> paras = new ArrayList<String>();//
		
		SuffixArraySequence sequence = new SuffixArraySequence(fasta, Constants.AMINO_ACIDS_18);
		SuffixArray sa = new SuffixArray(sequence);
		
		if(mode == 0){//HCD
			paras.add("/home/kwj/Dropbox/PAR/HCDT.par");
		}else if(mode == 1){//ETD
			paras.add("/home/kwj/Dropbox/PAR/ETDT.par");
		}else{
			paras.add("/home/kwj/Dropbox/PAR/HCDT.par");
			paras.add("/home/kwj/Dropbox/PAR/ETDT.par");
		}
		ArrayList<Tolerance> tols = new ArrayList<Tolerance>();//
		
		if(mode == 0){//HCD
			tols.add(new Tolerance(20f, true));
		}else if(mode == 1){//ETD
			tols.add(new Tolerance(0.5f, false));
		}else{
			tols.add(new Tolerance(20f, true));
			tols.add(new Tolerance(0.5f, false));
		}
		//
		Tolerance pmtol = new Tolerance(20f, true);
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		WindowFilter filter = new WindowFilter(6, 50);
		ArrayList<IPDGenerator> ipgs = new ArrayList<IPDGenerator>();
		

		for(int i=0; i<paras.size(); i++){
			ipgs.add(new IPDGenerator(paras.get(i), aaSet, i).filter(filter));
		}
		
		int cntr = 0, sn = 0;
		int matchedsn = 0, ptmMatchedsn = 0;
		float offsets[] = new float[300];
		
		for(File specFile : specFiles){
			if(!specFile.getName().endsWith(".mgf")) continue;
			
			Iterator<Spectrum> iterator = new SpectraIterator(specFile.getAbsolutePath(), new MgfSpectrumParser());
			ArrayList<Spectrum> paired = new ArrayList<Spectrum>();
			
			String title = "";
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				//if(spec.getCharge() != 2) continue;
				
				if(spec.getTitle().contains("experiment: 2"))
					title = spec.getTitle();
				
				if(mode == 0){//HCD
					if(spec.getTitle().contains("experiment: 2")){
						paired.add(spec);
					}else continue;
				}else if(mode == 1){//ETD
					if(spec.getTitle().contains("experiment: 1")){
						paired.add(spec);
					}else continue;
				}else{// HCD+ETD
					paired.add(spec);
					if(paired.size() < 2) continue;
				}
				
				sn++;
				SpectrumGraph graph = new SpectrumGraph(paired, ipgs, null, tols, pmtol, maxPTMNum > 0, true, false);
			
				UniNovo mg = new UniNovo(graph);
				
				ArrayList<DenovoReconstruction> out = new ArrayList<DenovoReconstruction>();
				
				float setAccuracy = mg.getOutput(out, numPerLength, minNumPerLength, minLength, maxLength, setAccuracyThreshold, maxGapNum, maxPTMNum, false);
				
				paired.clear();
				if(setAccuracy <= setAccuracyThreshold){
					System.out.println("# Filtered - " + specFile.getName() + "\t" + title + "\t" + setAccuracy);
					continue;
				}
				String s = ">>" + cntr++ + "\t";
				
				s += (specFile.getName() + "\t" + spec.getScanNum() + "\t" + title)+"\n";
				
				ArrayList<String> matches = new ArrayList<String>();
				
				for(DenovoReconstruction d : out){
					s += (d+"\t" + Math.round(d.getAccuracy()*100)+"%\t"+d.getLR())+"\n";
					for(String match : d.getMatches(sa, aaSet, false)){
						if(!matches.contains(match)) matches.add(match);
					}
				}
				
				
				if(matches.isEmpty()){
					ArrayList<Integer> off = new ArrayList<Integer>();
					for(DenovoReconstruction d : out){
						for(String match : d.getMatches(sa, aaSet, true)){
							int o = 50+Math.round(Float.parseFloat(match.split("\t")[1]));							
							//match = match.substring(0, match.indexOf('\t'));
							String[] t = match.split("\t");
							boolean isin = false;
							for(String prevmatch : matches){
								String[] pt = prevmatch.split("\t");
								if((pt[0]).equals(t[0])){
									isin = true;
									break;
								}
							}
							
							if(!isin){
								matches.add(match);
								off.add(o);
							}
						}
					}
					if(!matches.isEmpty()){
						ptmMatchedsn++;
						
						for(int o : off){
							offsets[o]++;
						}
					
					}
				}else matchedsn++;
				
				
				s += "Set Accuracy : " + Math.round(setAccuracy*100)+"%";
				
				if(!matches.isEmpty()){
					s += "\nMatches:";
					
					for(String match : matches){
						s += "\n" + match;
					}
				}
				System.out.println(s);
				op.println(s);
				
				//break;
			}
		//	if(ptmMatchedsn > 100) break;
		}
		String s = "\nTotal spec : " + sn + " Qualified spec : " + cntr + " Matched spec : " + matchedsn +
			" Matched spec with PTM : " + ptmMatchedsn + "\n";
		
		s += "offset = [\n";
		
		for(int i=0; i<offsets.length; i++){
			s += (i-50)+ "\t" + offsets[i] + "\n";
		}
		
		s+= "];";
		
		System.out.println(s);
		op.println(s);
		op.close();

	}

}
