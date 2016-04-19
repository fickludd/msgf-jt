package uninovoOld;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.BufferedLineReader;
import parser.MS2SpectrumParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import uninovoOld.DenovoReconstruction.DenovoReconstructionAccuracyComparator;
import uninovoOld.DenovoReconstruction.DenovoReconstructionLRComparator;





public class UniNovo {
	private SpectrumGraph graph;
	
	private void union(HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>> e, Node v, int l, int s, Edge edge){
		if(!e.containsKey(v)) e.put(v, new HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>());
		HashMap<Integer, HashMap<Integer, ArrayList<Edge>>> k = e.get(v);
		if(!k.containsKey(l)) k.put(l, new HashMap<Integer, ArrayList<Edge>>());
		HashMap<Integer, ArrayList<Edge>> q = k.get(l);
		if(!q.containsKey(s)) q.put(s, new ArrayList<Edge>());
		ArrayList<Edge> p = q.get(s);
		p.add(edge);
	}
	
	
	
	public UniNovo(SpectrumGraph graph){
		this.graph = graph;
	}
	
	
	public ArrayList<DenovoReconstruction> getCandidateDenovoReconstructions(int numPerLength, int minLength, int maxLength, int maxGapNum, int maxPTMNum, boolean accuracyFirst){
		ArrayList<DenovoReconstruction> candidateDenovoReconstructions = new ArrayList<DenovoReconstruction>();
		
		HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>> es = new HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>>();
		
		union(es, graph.getSourceNode(), 0, graph.getSourceNode().getLRScore(), new Edge(0));
		int minScore = Integer.MAX_VALUE;
		int maxScore = Integer.MIN_VALUE;
		for(Node v : graph.getNodes()){
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
			boolean onlyCGGap = false;
			
			if(maxGapNum < 0){
				maxGapNum = - maxGapNum;
				onlyCGGap = true;
				
			}
			for(int l=maxLength;l>=minLength;l--){
				ArrayList<DenovoReconstruction> candidates = new ArrayList<DenovoReconstruction>();
				for(int score = maxScore; score >= minScore; score--){
					backtrack(graph, candidateDenovoReconstructions, candidates, new ArrayList<Node> (), es, graph.getSinkNode(), l, score, maxGapNum, maxPTMNum, numPerLength, onlyCGGap, accuracyFirst);	
				}
			}
		}
		return candidateDenovoReconstructions;
	}
	
	
	
	private void backtrack(SpectrumGraph graph, ArrayList<DenovoReconstruction> candidateDenovoReconstructions,  ArrayList<DenovoReconstruction> rec, ArrayList<Node> suffix, HashMap<Node, HashMap<Integer, HashMap<Integer, ArrayList<Edge>>>> es, Node v, int l, int s,  int remainingGapNum, int remainingPTMNum, int numPerLength, boolean onlyCGGap, boolean accuracyFirst){
		if(rec.size() > numPerLength) return;
		HashMap<Integer, HashMap<Integer, ArrayList<Edge>>> subes = es.get(v);
		if(subes == null || subes.isEmpty()) return;
		
		HashMap<Integer, ArrayList<Edge>> ssubes = subes.get(l);
	
		if(ssubes == null || ssubes.isEmpty()) return;
		ArrayList<Edge> edges = ssubes.get(s);
		if(edges == null || edges.isEmpty()) return;
		
		if(v.equals(graph.getSourceNode())){// source
		
			suffix.add(v);
			DenovoReconstruction g = new DenovoReconstruction(suffix, graph);
			
			if(accuracyFirst){
				boolean isRedundant = false;
				for(int j=0; j< candidateDenovoReconstructions.size(); j++){
					DenovoReconstruction g1 = candidateDenovoReconstructions.get(j);
				
					if(g1.isRedundantWRT(g)){
						candidateDenovoReconstructions.remove(j--);
					}
				}
				
				for(int j=0; j< candidateDenovoReconstructions.size(); j++){
					DenovoReconstruction g1 = candidateDenovoReconstructions.get(j);
				
					if(g.isRedundantWRT(g1)){
						isRedundant = true;
						break;
					}
				}
				
				if(!isRedundant){
					rec.add(g);
					candidateDenovoReconstructions.add(g);
				}
				
			}else{		
				rec.add(g);
				if(!candidateDenovoReconstructions.contains(g))
					candidateDenovoReconstructions.add(g);
			}
			return;
		}

		ArrayList<Node> nextSuffix = new ArrayList<Node>();
		nextSuffix.addAll(suffix);
		nextSuffix.add(v);
		
		for(Edge edge : edges){
			if(onlyCGGap){
				if(edge.getLeftNode().getMass() != 0 && !edge.getRightNode().isSink() && edge.isGap()) continue;
			}
			
			int nextRemainingGapNum = remainingGapNum;
			if(edge.isGap()) nextRemainingGapNum--;
			if(nextRemainingGapNum<0) continue;
			
			int nextRemainingPTMNum = remainingPTMNum;
			if(edge.isPTM()) nextRemainingPTMNum--;
			if(nextRemainingPTMNum<0) continue;

			backtrack(graph, candidateDenovoReconstructions, rec, nextSuffix, es, edge.getLeftNode(), l!=0? l-1 : l, s-edge.getLRScore(), nextRemainingGapNum, nextRemainingPTMNum, numPerLength, onlyCGGap, accuracyFirst);
		}
	}
	
	
	
	
	public float selectOutputFromCandidates(ArrayList<DenovoReconstruction> out, ArrayList<DenovoReconstruction> candidateDenovoReconstructions, int numPerLength, int minNumPerLength, float setAccuracyThreshold, int maxPTMNum, boolean accuracyFirst){
		float acc = 0;
		boolean met = false;
		if(!accuracyFirst){
			Collections.sort(candidateDenovoReconstructions, Collections.reverseOrder(DenovoReconstructionLRComparator.get()));

			if(setAccuracyThreshold <=0){
				for(int i=0;i<Math.min(candidateDenovoReconstructions.size(), numPerLength);i++)
					out.add(candidateDenovoReconstructions.get(i));
			}else{			
				int maxGn = 0;
				int maxL = 0;
				for(int i=0; i<candidateDenovoReconstructions.size() ; i++){
					maxL = Math.max(maxL, candidateDenovoReconstructions.get(i).length());
					maxGn = Math.max(maxGn, candidateDenovoReconstructions.get(i).getGapNum());
				}
				
				for(int l = maxL;l>=0;l--){
					float sumAcc = 0;
					for(int gn = 0; gn <= maxGn; gn++){
						out.clear();
						sumAcc = 0;
						for(int i=0; i<candidateDenovoReconstructions.size() ; i++){
							DenovoReconstruction g = candidateDenovoReconstructions.get(i);
							if(g.length() <= l && g.getGapNum() <= gn){
								out.add(g);
								sumAcc += g.getAccuracy();
							}else continue;
							if(out.size() >= Math.min(candidateDenovoReconstructions.size(), numPerLength))break;
						}
						if(sumAcc < setAccuracyThreshold) continue;
						
						acc = getSetAccuracy(out, maxPTMNum);
						
						if(out.size() >= Math.min(candidateDenovoReconstructions.size(), minNumPerLength) && acc >= setAccuracyThreshold){
							met = true;
							break;
						}
					}
					if(met) break;
				}
			}
		}else{
			Collections.sort(candidateDenovoReconstructions, Collections.reverseOrder(DenovoReconstructionAccuracyComparator.get()));
			for(int i=0; i<Math.min(candidateDenovoReconstructions.size(), numPerLength) ; i++){
				DenovoReconstruction g = candidateDenovoReconstructions.get(i);

				out.add(g);
				acc = getSetAccuracy(out, maxPTMNum);
				if(out.size() >= Math.min(candidateDenovoReconstructions.size(), minNumPerLength) && acc >= setAccuracyThreshold) break;
			}			
		}
		return Math.max(0f, acc);
	}
	
	public float getOutput(ArrayList<DenovoReconstruction> out, int numPerLength, int minNumPerLength, int minLength, int maxLength, float setAccuracyThreshold, int maxGapNum, int maxPTMNum, boolean accuracyFirst){
		float setAccuracy = selectOutputFromCandidates(out, getCandidateDenovoReconstructions(Math.max(100, numPerLength * 3), minLength, maxLength, maxGapNum, maxPTMNum, accuracyFirst), numPerLength, minNumPerLength, setAccuracyThreshold, maxPTMNum, accuracyFirst);
		return setAccuracy;
	}
	
	private float getSetAccuracy(ArrayList<DenovoReconstruction> recs, int maxPTMNum){

		ArrayList<DenovoReconstruction> toConsider = new ArrayList<DenovoReconstruction>();
		
		new ArrayList<DenovoReconstruction>();
		for(int k=0; k< recs.size(); k++){
			DenovoReconstruction g1 = recs.get(k);
			boolean isRedundant = false;
			for(int l=0; l< recs.size(); l++){
				if(k == l) continue;
				DenovoReconstruction g2 = recs.get(l);
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
			
			DenovoReconstruction g = toConsider.get(k);
			DenovoReconstruction g2 = toConsider.get(k-1);
			
			DenovoReconstruction u = g.getUnionDenovoReconstruction(g2, maxPTMNum);
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
	
	static private void printUsageAndExit(){
		System.out.println("Usage : java -jar -Xmx2000m UniNovo.jar" +
				"\n -i input mgf,mzXML, or ms2 files seperated by :" +
				"\n -o output file prefix" +
				"\n -t ion tolerances seperated by : (each ending in ppm or Da)" +
				"\n -pt precursor ion tolerance (ending in ppm or Da)" +
				"\n -f fragmentation methods seperated by : (CID/ETD/HCD)" +
				
				//"\n -MSGFDB MSGFDB output file name. If specified, the identified spectra will be skipped (with PepFDR > 0.01)" +
				"\n -e enzyme applied (0: No enzyme specificity, 1: Trypsin (default), 2: LysC, 3: AspN)" +
				"\n -c number of 13C (default : 0)" +
				"\n -l minimum length of reconstructions (default : 5)");
			
		System.out.println(" -acc set accuracy threshold (0.0-0.9 : default 0.8)" +
				"\n -n number of de novo sequences per one spectrum (1-200 : default 100)" +
				"\n -g number of possible mass gaps per each sequence (2-10 : default 10)" +
				//"\n -Nmod number of possible blind modifications per each peptide (0-3 : default 0)" +
				"\n -par use user trained parameter file");
		
		System.out.println("Example : " +
				"\n  java -jar -Xmx2000m UniNovo.jar -i /home/test.mgf -o /home/testOutput -t 0.1Da -pt 20ppm -f CID -e 1" +
				"\n java -jar -Xmx2000m UniNovo.jar -i CID.mgf:ETD.mgf -o CIDETDOutput -t 0.1Da:0.1Da -pt 20ppm -f CID:ETD" +
				"\n java -jar -Xmx2000m UniNovo.jar -i CID.mgf:ETD.mgf:HCD.mgf -o CIDETDHCDOutput -t 0.1Da:0.1Da:20ppm -pt 20ppm -f CID:ETD:HCD");
		
		System.exit(0);
	}
	
	public static void main(String[] args) throws IOException{
		
		if(args.length < 1) printUsageAndExit();
		
		ArrayList<String> specFiles = new ArrayList<String>();//
		ArrayList<String> paras = new ArrayList<String>();//
		ArrayList<Tolerance> tols = new ArrayList<Tolerance>();//
		HashSet<String> identifiedByMSGFDB = new HashSet<String>();
		
		Tolerance pmtol = null;//
		Tolerance inputPmtol = null;
		
		String prefix = "";
		String MSGFDBoutfile = "";
		Enzyme enzyme = Enzyme.TRYPSIN;
		
		String fasta = null;
		
		int numRecs = 100;
		int minNumRecs = 0;
		int minLength = 5;
		int maxLength = 60;
		float setAccuracyThreshold = 0.800f;
		int maxGapNum = 10;
		int maxPTMNum = 0;
		int num13C = 1;
		int searchMod = 0; // 0 exact 1 blind 2 mut 
		float FDR = 0.01f;
		boolean accuracyFirst = false;
		//boolean reuse = false;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		int maxSN = Integer.MAX_VALUE;
		
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
			else if(arg.startsWith("-Nmod")) mode = 9;
			else if(arg.startsWith("-acc")) mode = 10;
			else if(arg.startsWith("-aa")) mode = 11;
			else if(arg.startsWith("-c")) mode = 12;
			else if(arg.startsWith("-MSGFDB")) mode = 13;
			else if(arg.startsWith("-s")) mode = 14;
			else if(arg.startsWith("-l")) mode = 15;
		//	else if(arg.startsWith("-u")) reuse = true;
			else if(arg.startsWith("-par")) mode = 16;
			else if(arg.startsWith("-maxSN")) mode = 17; // hidden
			else{
				if(mode < 0) continue;
				if(mode == 0){
					for(String file : arg.split(":"))
						specFiles.add(file);
				}else if(mode == 1){
					prefix = arg;
				}else if(mode == 2){
					for(String tol : arg.split(":")){
						tols.add(Tolerance.parseToleranceStr(tol));
					}
				}else if(mode == 3){
					inputPmtol = Tolerance.parseToleranceStr(arg);
					if(inputPmtol.getToleranceAsDa(500f) < 0.5f)
						pmtol = inputPmtol;
					else 
						pmtol = new Tolerance(0.5f);
				}else if(mode == 4){
					for(String parafile : arg.split(":")){
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
				}else if(mode == 6){
					fasta = arg;
				}else if(mode == 7){
					numRecs = Integer.parseInt(arg);
				}else if(mode == 8){
					maxGapNum = Integer.parseInt(arg);
				}else if(mode == 9){
					maxPTMNum = Integer.parseInt(arg);
				}else if(mode == 10){
					setAccuracyThreshold = Float.parseFloat(arg);
				}else if(mode == 11){
					aaSet = AminoAcidSet.getAminoAcidSetFromModFile(arg);
				}else if(mode == 12){
					num13C = Integer.parseInt(arg);
					SpectrumGraph.setNum13C(num13C);
				}else if(mode == 13){
					MSGFDBoutfile = arg;
				}else if(mode == 14){
					searchMod = Integer.parseInt(arg);		
					if(searchMod == 1) maxPTMNum = 1;
				}else if(mode == 15){
					minLength = Integer.parseInt(arg);		
				}else if(mode == 16){
					for(String par : arg.split(":")){
						paras.add(par);
					}
				}else if(mode == 17){
					maxSN = Integer.parseInt(arg);		
				}
			}
		}
		
		String deNovoOutput = prefix + ".den";
		String grcOutput = prefix +".grc";
		String scoringOutput = prefix + ".mic";
		
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
		
		ArrayList<IPDGenerator> ipgs = new ArrayList<IPDGenerator>();
		
		WindowFilter filter = new WindowFilter(6, 50);
	
		for(int i=0; i<specFiles.size(); i++){
			ipgs.add(new IPDGenerator(paras.get(i), aaSet, i).filter(filter));
		}
		
		if(fasta != null)
			accuracyFirst = true;
		
		if(!MSGFDBoutfile.isEmpty()){
			String s;
			BufferedLineReader in = new BufferedLineReader(MSGFDBoutfile);
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");
				if(Float.parseFloat(token[14]) > FDR) continue;
				
				identifiedByMSGFDB.add(token[1]+"\t"+token[4]);
			}
			in.close();
		}
		
		
		ArrayList<Iterator<Spectrum>> iterators = new ArrayList<Iterator<Spectrum>>();
			
		for(String file : specFiles){
			iterators.add(getSpectralIterator(file));
		}
	
		PrintStream deNovoFile = new PrintStream(deNovoOutput);
		PrintStream grcFile = null;
		PrintStream miscFile = null;
		
		if(fasta != null){
			grcFile = new PrintStream(grcOutput);
			miscFile = new PrintStream(scoringOutput);
		}
		
		int cumCount = 0;
		int qsn = 0;
		int correctpm = 0;
		
		boolean first = true;
		
		if(fasta != null){
			for(Tolerance tol : tols){
				miscFile.println(">>TOL\t" + tol);
			}
			miscFile.println(">>PMTOL\t" + pmtol);
		}
		
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
				
				if(identifiedByMSGFDB.contains((spec.getScanNum()+1)+"\t"+spec.getPrecursorPeak().getMz())){
					isIdentified = true; break;
				}
				
				if(spec.getCharge() == 0) spec.setCharge(2);
				//spec.correctParentMass();
				spectra.add(spec);
			}
			
			if(isIdentified || spectra.isEmpty()) continue;

			SpectrumGraph graph = new SpectrumGraph(spectra, ipgs, enzyme, tols, inputPmtol, maxPTMNum > 0, true, true);
					
			if(graph.getSinkNode().isCorrect(spectra.get(0).getAnnotation())){
				correctpm++;
			}
			UniNovo adaNovo = new UniNovo(graph);
			
			ArrayList<DenovoReconstruction> out = new ArrayList<DenovoReconstruction>();
			float setAccuracy = adaNovo.getOutput(out, numRecs, minNumRecs, minLength, maxLength, setAccuracyThreshold, maxGapNum, maxPTMNum, accuracyFirst);
			
			if(fasta != null){
				if(first)
					first = false;
				else
					miscFile.println("***");
			
				String miscHeaderLine = String.format("#%d\t%d\t%s\t%.5f\t", cumCount, spectra.get(0).getScanNum(), specFiles.get(0), spectra.get(0).getParentMass());
				
				for(Spectrum s : spectra){
					miscHeaderLine += s.getCharge()+",";
				}
				miscFile.println(miscHeaderLine);
				miscFile.println(graph.toFileString());
			}
			
			System.out.print(">> " + spectra.get(0).getScanNum() + "\t" +pm + "\t" + origCharge + "\t" + spectra.get(0).getActivationMethod() +"\t"+ setAccuracy);

			if(setAccuracy < setAccuracyThreshold || out.isEmpty()){
				System.out.println("\tfiltered out " + out.size());
				continue;
			}else System.out.println();
			deNovoFile.println(">>\t" + spectra.get(0).getScanNum() + "\t" + spectra.get(0).getTitle() + "\t" +pm + "\t" + origCharge);
			
			if(fasta != null){
				String grcHeaderLine = String.format("#%d\t%d\t%s\t%.5f\t%d", cumCount, spectra.get(0).getScanNum(), specFiles.get(0), pm, origCharge);
				grcFile.println(grcHeaderLine);
			}
			
			deNovoFile.println("Reconstruction\tScore\tAccuracy");
			
			for(DenovoReconstruction r : out){
				Peptide cp = spectra.get(0).getAnnotation();
				
				deNovoFile.print(r + "\t" + String.format("%.2f", r.getLR()) + "\t" + String.format("%.2f %%", r.getAccuracy()*100));
				if(cp!=null && r.isCorrect(cp)){
					deNovoFile.print("\t*");
				}
				deNovoFile.println();
				if(fasta != null){
					for(ArrayList<Integer> gp : r.getGapMassRepresentation(aaSet))
						grcFile.println(gp);					
				}
			}
			
			deNovoFile.println("Set Accuracy: " + setAccuracy);
			qsn++;
			if(qsn > maxSN) break;
		}
		if(fasta != null) {
			miscFile.println("*#*");
			grcFile.close();
			miscFile.close();
		}
		
		deNovoFile.close();
		
		System.out.println("Qualified spectra: " + qsn + "\nTotal spectra: " + cumCount + "\t" + correctpm);
		
	}
	
}
