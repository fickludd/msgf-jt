/***************************************************************************
  * Title:          
  * Author:         Kyowon Jeong
  * Last modified:  
  *
  * Copyright (c) 2008-2009 The Regents of the University of California
  * All Rights Reserved
  * See file LICENSE for details.
  ***************************************************************************/
package uninovo;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import uninovo.DenovoSequence.DenovoSequenceLRComparator;
import uninovo.parser.MS2SpectrumParser;
import uninovo.parser.MgfSpectrumParser;
import uninovo.parser.MzXMLSpectraIterator;
import uninovo.util.AminoAcidSet;
import uninovo.util.Enzyme;
import uninovo.util.SpectraIterator;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
import uninovo.util.WindowFilter;


/**
 * UniNovo class takes spectrum files and outputs de novo seqeuences. 
 * @author kyowon
 */

public class UniNovo {

	/**
	 * Backtrack the dynamic programming table to get sequences with a fixed length
	 * @param graph spectrum graph
	 * @param candidateSequences output candidate sequences
	 * @param seqs sequences with a fixed length  
	 * @param suffix suffix of the sequences
	 * @param table dynamic programming table
	 * @param v current node
	 * @param l remaining length of sequences
	 * @param s score of the current edge
	 * @param allowedGapNum number of allowed gaps
	 * @param numPerLength allowed number of sequences per length
	 */
	static private void backtrack(SpectrumGraph graph, ArrayList<DenovoSequence> candidateSequences,  ArrayList<DenovoSequence> seqs, ArrayList<Node> suffix, HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>> table, Node v, int l, int s,  int allowedGapNum, int numPerLength){
		if(seqs.size() > numPerLength) return;
		HashMap<Integer, HashMap<Integer, ArrayList<Edge>>> subes = table.get(v);
		if(subes == null || subes.isEmpty()) return;
		
		HashMap<Integer, ArrayList<Edge>> ssubes = subes.get(l);
	
		if(ssubes == null || ssubes.isEmpty()) return;
		ArrayList<Edge> edges = ssubes.get(s);
		if(edges == null || edges.isEmpty()) return;
		
		if(v.equals(graph.getSourceNode())){// source
		
			suffix.add(v);
			DenovoSequence g = new DenovoSequence(suffix, graph);
			
			boolean isRedundant = false;
			for(int j=0; j< candidateSequences.size(); j++){
				DenovoSequence g1 = candidateSequences.get(j);
			
				if(g1.isRedundantWRT(g)){
					candidateSequences.remove(j--);
				}
			}
			
			for(int j=0; j< candidateSequences.size(); j++){
				DenovoSequence g1 = candidateSequences.get(j);
			
				if(g.isRedundantWRT(g1)){
					isRedundant = true;
					break;
				}
			}
			
			if(!isRedundant){
				seqs.add(g);
				candidateSequences.add(g);
			}
		
			return;
		}

		ArrayList<Node> nextSuffix = new ArrayList<Node>();
		nextSuffix.addAll(suffix);
		nextSuffix.add(v);
		
		for(Edge edge : edges){						
			int nextRemainingGapNum = allowedGapNum;
			if(edge.isGap()) nextRemainingGapNum--;
			if(nextRemainingGapNum<0) continue;			

			backtrack(graph, candidateSequences, seqs, nextSuffix, table, edge.getLeftNode(), l!=0? l-1 : l, s-edge.getLRScore(), nextRemainingGapNum, numPerLength);
		}
	}
	
	/**
	 * Generate de novo sequences 
	 * @param graph spectrum graph
	 * @param out the output sequences
	 * @param numSeq number of sequences 
	 * @param minNumPerLength min number of sequences per length 
	 * @param minLength minimum length of sequences
	 * @param maxLength maximum length of sequences
	 * @param setAccuracyThreshold set accuracy threshold
	 * @param maxGapNum maximum number of gaps allowed
	 * @return set accuracy
	 */
	static public float generateDenovoSequences(SpectrumGraph graph, ArrayList<DenovoSequence> out, int numSeq, int minNumPerLength, int minLength, int maxLength, float setAccuracyThreshold, int maxGapNum){
		float setAccuracy = selectFromCandidateSequences(out, getCandidateDenovoSequences(graph, Math.max(100, numSeq * 3), minLength, maxLength, maxGapNum), numSeq, minNumPerLength, setAccuracyThreshold);
		return setAccuracy;
	}
	
	/**
	 * Generate candidate sequences on spectrum graph using a simple dynamic programming
	 * @param graph spectrum graph
	 * @param numPerLength the number of sequences generated per each length
	 * @param minLength sequences shorter than minLength will not be generated
	 * @param maxLength sequences longer than maxLength will not be generated
	 * @param maxGapNum sequences with more gaps then maxGapNum will not be generated
	 * @return the candidate sequences
	 */
	static public ArrayList<DenovoSequence> getCandidateDenovoSequences(SpectrumGraph graph, int numPerLength, int minLength, int maxLength, int maxGapNum){
		ArrayList<DenovoSequence> candidateSequences = new ArrayList<DenovoSequence>();
		HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>> es = new HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>>();
		
		union(es, graph.getSourceNode(), 0, graph.getSourceNode().getLRScore(), new Edge(0));
		int minScore = Integer.MAX_VALUE;
		int maxScore = Integer.MIN_VALUE;
		for(Node v : graph.getFilteredNodes()){
			for(Edge e: graph.getLinkedEdges(v)){
				Node ln = e.getLeftNode();
				HashMap<Integer, HashMap<Integer, ArrayList<Edge>>> subes = es.get(ln);
				if(subes == null || subes.isEmpty()) continue;
				for(int l : subes.keySet()){
					HashMap<Integer, ArrayList<Edge>> ssubes = subes.get(l);
					for(int s : ssubes.keySet()){
						int ns = s+e.getLRScore();
						union(es, v, l+1, ns, e);
						if(v.equals(graph.getSinkNode())){
							minScore = Math.min(minScore, ns);
							maxScore = Math.max(maxScore, ns);
						}
					}
				}
			}
		
			for(int l=maxLength;l>=minLength;l--){
				ArrayList<DenovoSequence> rec = new ArrayList<DenovoSequence>();
				for(int score = maxScore; score >= minScore; score--){
					backtrack(graph, candidateSequences, rec, new ArrayList<Node> (), es, graph.getSinkNode(), l, score, maxGapNum, numPerLength);	
				}
			}
		}
		
		return candidateSequences;
	}
	
	/**
	 * Calculate the set accuracy of a set of sequences
	 * @param seqs the set of sequences
	 * @return the set accuracy
	 */
	static private float getSetAccuracy(ArrayList<DenovoSequence> seqs){

		ArrayList<DenovoSequence> toConsider = new ArrayList<DenovoSequence>();
		
		new ArrayList<DenovoSequence>();
		for(int k=0; k< seqs.size(); k++){
			DenovoSequence g1 = seqs.get(k);
			boolean isRedundant = false;
			for(int l=0; l< seqs.size(); l++){
				if(k == l) continue;
				DenovoSequence g2 = seqs.get(l);
				if(g1.isRedundantWRT(g2)){
					isRedundant = true;
					break;
				}
			}
			if(!isRedundant)
				toConsider.add(g1);
		}

		int i=toConsider.size()-1;
		if(i<0) return 0;
	
		float p=1-toConsider.get(0).getAccuracy();
		float sumProb = 1-p;
		float acc = sumProb;
		
		for(int k=1;k<=i;k++){
			
			DenovoSequence g = toConsider.get(k);
			DenovoSequence g2 = toConsider.get(k-1);
			
			DenovoSequence u = g.getUnionDenovoSequence(g2);
			float up = g.getAccuracy() + g2.getAccuracy() - u.getAccuracy();
			up = Math.max(up, 0.85f*g.getAccuracy());
			up = Math.max(up, 0.85f*g2.getAccuracy());
			
			float f = (1 - up)/(1-g2.getAccuracy());
			f= Math.max(f, 0);
			f = Math.min(1f, f);					
			p *= f;			
			sumProb += g.getAccuracy();
			acc = Math.max(acc, g.getAccuracy());
			
		}
		
		acc = Math.max(acc, Math.min(sumProb, 1-p));
		acc = Math.min(acc, 1);
		return acc;
	}
	
	/**
	 * Get an iterator of spectra 
	 * @param spectrumFileName spectrum file name
	 * @return iterator of spectra
	 */
	static public Iterator<Spectrum> getSpectralIterator(String spectrumFileName){
		Iterator<Spectrum> iterator = null;
		try {
			if(spectrumFileName.endsWith("mgf"))
				iterator = new SpectraIterator(spectrumFileName, new MgfSpectrumParser());
			else if (spectrumFileName.endsWith("mzXML")) 
				iterator = new MzXMLSpectraIterator(spectrumFileName);
			else if (spectrumFileName.endsWith("ms2")) 
			  iterator = new SpectraIterator(spectrumFileName, new MS2SpectrumParser());
		}catch (IOException e) {
			 System.err.println("IOException no spectrum file found named " + spectrumFileName);
			 e.printStackTrace();
			 System.exit(-1);
		}
		return iterator;
	}
	
	/**
	 * The main function for UniNovo.
	 *
	 * @param args specifies the input parameters
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void main(String[] args) throws IOException{
		
		if(args.length < 1) printUsageAndExit();
		
		ArrayList<String> specFiles = new ArrayList<String>();//
		ArrayList<String> paras = new ArrayList<String>();//
		ArrayList<Tolerance> tols = new ArrayList<Tolerance>();//
		Tolerance inputPmtol = null;
		
		String prefix = "";
		Enzyme enzyme = Enzyme.TRYPSIN;
		
		int numRecs = 100;
		int minNumRecs = 0;
		int minLength = 5;
		int maxLength = 60;
		float setAccuracyThreshold = 0.800f;
		int maxGapNum = 10;
		int num13C = 0;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
			
		ArrayList<String> frags = new ArrayList<String>();
		int mode = -1;
		for(String arg : args){
			if(arg.startsWith("-i")) mode = 0;
			else if(arg.startsWith("-o")) mode = 1;
			else if(arg.startsWith("-t")) mode = 2;
			else if(arg.startsWith("-pt")) mode = 3;
			else if(arg.startsWith("-f")) mode = 4;
			else if(arg.startsWith("-e")) mode = 5;
			else if(arg.startsWith("-d")) mode = 6;
			else if(arg.startsWith("-n")) mode = 7;
			else if(arg.startsWith("-g")) mode = 8;
			else if(arg.startsWith("-acc")) mode = 10;
			else if(arg.startsWith("-c")) mode = 12;
			else if(arg.startsWith("-l")) mode = 15;
			else if(arg.startsWith("-par")) mode = 16;
			else{
				if(mode < 0) continue;
				if(mode == 0){
					for(String file : arg.split(","))
						specFiles.add(file);
				}else if(mode == 1){
					prefix = arg;
				}else if(mode == 2){
					for(String tol : arg.split(",")){
						tols.add(Tolerance.parseToleranceStr(tol));
					}
				}else if(mode == 3){
					inputPmtol = Tolerance.parseToleranceStr(arg);					
				}else if(mode == 4){
					for(String parafile : arg.split(",")){
						frags.add(parafile);
					}
				}else if(mode == 5){
					int e = Integer.parseInt(arg);
					if(e==0){
						enzyme = null;
					}else if(e==1){
						enzyme = Enzyme.TRYPSIN;
					}else if(e==2){
						enzyme = Enzyme.LysC;
					}else if(e==3)
						enzyme = Enzyme.AspN;			
				}else if(mode == 7){
					numRecs = Integer.parseInt(arg);
				}else if(mode == 8){
					maxGapNum = Integer.parseInt(arg);
				}else if(mode == 10){
					setAccuracyThreshold = Float.parseFloat(arg);
				}else if(mode == 12){
					num13C = Integer.parseInt(arg);
					SpectrumGraph.setNum13C(num13C);
				}else if(mode == 15){
					minLength = Integer.parseInt(arg);		
				}else if(mode == 16){
					for(String par : arg.split(",")){
						paras.add(par);
					}
				}
			}
		}
		
		String deNovoOutput = prefix + ".den";
	
		if(paras.isEmpty()){
			for(int i=0;i<frags.size();i++){
				String frag = frags.get(i);
				String suffix = "_train.pars";
				
				if(!frag.equals("HCD")){
					if(enzyme != null && enzyme.equals(Enzyme.LysC)){
						suffix = "LysC"+suffix;
					}if(enzyme != null && enzyme.equals(Enzyme.AspN)){
						suffix = "AspN"+suffix;
					}else{
						suffix = "Trypsin"+suffix;
					}
					
					if(tols.get(i).getToleranceAsDa(500f)<0.02f) suffix = "FT"+suffix;
				}
				paras.add("Pars" + System.getProperty("file.separator") +frag+suffix);
			}
		}
		
		ArrayList<PeakIonProbabilityGenerator> ipgs = new ArrayList<PeakIonProbabilityGenerator>();
		
		WindowFilter filter = new WindowFilter(6, 50);
	
		for(int i=0; i<specFiles.size(); i++){
			ipgs.add(new PeakIonProbabilityGenerator(paras.get(i), aaSet, i).setFilter(filter));
		}
		
		ArrayList<Iterator<Spectrum>> iterators = new ArrayList<Iterator<Spectrum>>();
			
		for(String file : specFiles){
			iterators.add(getSpectralIterator(file));
		}
	
		PrintStream deNovoFile = new PrintStream(deNovoOutput);
		int cumCount = 0;
		int qsn = 0;		
		
		while(iterators.get(0).hasNext()){
			ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
			
			float pm = 0;
			int origCharge = 0;
			
			cumCount++;
			
			boolean isIdentified = false;
			
			for(int i=0; i<iterators.size(); i++){
				Iterator<Spectrum> iterator = iterators.get(i);
				Spectrum spec = iterator.next();

				if(pm == 0){
					pm = spec.getPrecursorPeak().getMz();
					origCharge = spec.getCharge();
				}
				else if(pm != spec.getPrecursorPeak().getMz()) continue;
				if(spec.getCharge() == 0) spec.setCharge(2);	
				//spec.correctParentMass();
				spectra.add(spec);
			}	
			if(isIdentified || spectra.isEmpty()) continue;
			SpectrumGraph graph = new SpectrumGraph(spectra, ipgs, enzyme, tols, inputPmtol, false, true, true);		
			
			ArrayList<DenovoSequence> out = new ArrayList<DenovoSequence>();
			float setAccuracy = UniNovo.generateDenovoSequences(graph, out, numRecs, minNumRecs, minLength, maxLength, setAccuracyThreshold, maxGapNum);
			
			System.out.print(">> " + spectra.get(0).getScanNum() + "\t" +pm + "\t" + origCharge + "\t" + spectra.get(0).getActivationMethod() +"\t"+ setAccuracy);

			if(setAccuracy < setAccuracyThreshold || out.isEmpty()){
				System.out.println("\tfiltered out " + out.size());
				continue;
			}else System.out.println();
			deNovoFile.println(">>\t" + spectra.get(0).getScanNum() + "\t" + spectra.get(0).getTitle() + "\t" +pm + "\t" + origCharge);			
			deNovoFile.println("Reconstruction\tScore\tAccuracy");
			
			for(DenovoSequence r : out){
				deNovoFile.print(r + "\t" + String.format("%.2f", r.getLR()) + "\t" + String.format("%.2f %%", r.getAccuracy()*100));
				deNovoFile.println();				
			}
			
			deNovoFile.println("Set Accuracy: " + setAccuracy);
			qsn++;
		
		}
		deNovoFile.close();
		
		System.out.println("Qualified spectra: " + qsn + "\nTotal spectra: " + cumCount);
		
	}
	
	/**
	 * Print usage and exit
	 */
	static private void printUsageAndExit(){
		System.out.println("Usage : java -jar -Xmx2000m UniNovo.jar" +
				"\n -i input mgf,mzXML, or ms2 files seperated by comma" +
				"\n -o output file prefix" +
				"\n -t ion tolerances seperated by comma (each ending in ppm or Da)" +
				"\n -pt precursor ion tolerance (ending in ppm or Da)" +
				"\n -f fragmentation methods seperated by comma (CID/ETD/HCD)" +
				"\n -e enzyme applied (0: No enzyme specificity, 1: Trypsin (default), 2: LysC, 3: AspN)" +
				"\n -c number of 13C (default : 0)" +
				"\n -l minimum length of reconstructions (default : 5)");
			
		System.out.println(" -acc set accuracy threshold (0.0-0.9 : default 0.8)" +
				"\n -n number of de novo sequences per one spectrum (1-200 : default 100)" +
				"\n -g number of possible mass gaps per each sequence (2-10 : default 10)" +
				"\n -par use user trained parameter file(s), seperated by comma");
		
		System.out.println("Example : " +
				"\n  java -jar -Xmx2000m UniNovo.jar -i /home/test.mgf -o /home/testOutput -t 0.1Da -pt 20ppm -f CID -e 1" +
				"\n java -jar -Xmx2000m UniNovo.jar -i CID.mgf,ETD.mgf -o CIDETDOutput -t 0.1Da,0.1Da -pt 20ppm -f CID,ETD" +
				"\n java -jar -Xmx2000m UniNovo.jar -i CID.mgf,ETD.mgf,HCD.mgf -o CIDETDHCDOutput -t 0.1Da,0.1Da,20ppm -pt 20ppm -f CID,ETD,HCD");
		
		System.exit(0);
	}
	
	/**
	 * Select sequences to meet set accuracy threshold from candidate sequences
	 * @param selectedSequences the selected sequences
	 * @param candidateSequences candidate sequences
	 * @param numSeq number of sequences
	 * @param minNumPerLength min number of sequences per length 
	 * @param setAccuracyThreshold set accuracy threshold
	 * @return the set accuracy of the set of selected sequences
	 */
	static public float selectFromCandidateSequences(ArrayList<DenovoSequence> selectedSequences, ArrayList<DenovoSequence> candidateSequences, int numSeq, int minNumPerLength, float setAccuracyThreshold){
		float acc = 0;
		boolean met = false;
		Collections.sort(candidateSequences, Collections.reverseOrder(DenovoSequenceLRComparator.get()));

		if(setAccuracyThreshold <=0){
			for(int i=0;i<Math.min(candidateSequences.size(), numSeq);i++)
				selectedSequences.add(candidateSequences.get(i));
			acc = getSetAccuracy(selectedSequences);
		}else{			
			int maxGn = 0;
			int maxL = 0;
			for(int i=0; i<candidateSequences.size() ; i++){
				maxL = Math.max(maxL, candidateSequences.get(i).length());
				maxGn = Math.max(maxGn, candidateSequences.get(i).getGapNum());
			}
			
			for(int l = maxL;l>=0;l--){
				float sumAcc = 0;
				for(int gn = 0; gn <= maxGn; gn++){
					selectedSequences.clear();
					sumAcc = 0;
					for(int i=0; i<candidateSequences.size() ; i++){
						DenovoSequence g = candidateSequences.get(i);
						if(g.length() <= l && g.getGapNum() <= gn){
							selectedSequences.add(g);
							sumAcc += g.getAccuracy();
						}else continue;
						if(selectedSequences.size() >= Math.min(candidateSequences.size(), numSeq))break;
					}
					if(sumAcc < setAccuracyThreshold) continue;
					
					acc = getSetAccuracy(selectedSequences);
					
					if(selectedSequences.size() >= Math.min(candidateSequences.size(), minNumPerLength) && acc >= setAccuracyThreshold){
						met = true;
						break;
					}
				}
				if(met) break;
			}
		}
		return Math.max(0f, acc);
	}
	
	/**
	 * A subfunction for dynamic programming to generate sequences on the spectrum graph
	 * @param table dynamic programming table
	 * @param v current node
	 * @param l remaining length of sequences
	 * @param s score of the added edge
	 * @param edge edge that will be added to table
	 */
	static private void union(HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>> table, Node v, int l, int s, Edge edge){
		if(!table.containsKey(v)) table.put(v, new HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>());
		HashMap<Integer, HashMap<Integer, ArrayList<Edge>>> k = table.get(v);
		if(!k.containsKey(l)) k.put(l, new HashMap<Integer, ArrayList<Edge>>());
		HashMap<Integer, ArrayList<Edge>> q = k.get(l);
		if(!q.containsKey(s)) q.put(s, new ArrayList<Edge>());
		ArrayList<Edge> p = q.get(s);
		p.add(edge);
	}
}
