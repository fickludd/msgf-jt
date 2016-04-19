package uninovoOld.analysis;

import java.io.File;
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
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;
import uninovoOld.DenovoReconstruction;
import uninovoOld.Edge;
import uninovoOld.IPDGenerator;
import uninovoOld.Node;
import uninovoOld.SpectrumGraph;
import uninovoOld.UniNovo;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.train.PeakGenerator;

public class AccuracyTest {
		
	public static void main(String[] args) throws IOException {
		
	String[] frag = {
			//"HCD",
			"CID"
		};
		
		String enz = "Trypsin";
		
		int charge = 2;

		int numPerLength = 20;

		int maxPTMNum = 0;
		
		float setAccuracyThreshold = 0.0f;
	
		if(args.length > 1){
			frag = new String[1];
			frag[0] = args[0];
			enz = args[1];
			
			charge = Integer.parseInt(args[2]);
			numPerLength = Integer.parseInt(args[3]);
			setAccuracyThreshold = Float.parseFloat(args[4]);
		}
		
		Enzyme enzyme = Enzyme.TRYPSIN;
		
		String inputmgf[] = {
			//	"/home/kwj/workspace/inputs/Training/MSGFDB_Tryp_FT_CID.mgf",
			System.getProperty("user.home") +  "/Dropbox/Test/"+frag[0] + (frag[0].equals("HCD")? "" : enz) +"_test.mgf"
			//"/home/kwj/Dropbox/Test/Standard1388_corrected.mgf"
				//"/home/kwj/Dropbox/Test/Modified/test.mgf"
				//	"/home/kwj/workspace/inputs/Standard1388.mgf"
				//"/home/kwj/workspace/inputs/MultiEnzymes/CID_IT_AspN/MSGFDB_AspN_F3.mgf",
		//"/home/kwj/Dropbox/Test/HCDETD/HCD.mgf",
	//	"/home/kwj/Dropbox/Test/HCDETD/ETD.mgf"
		//	"/home/kwj/Dropbox/Test/CIDETD/CID.mgf",
		//	"/home/kwj/Dropbox/Test/CIDETD/ETD.mgf"
		};
		
		int maxGapNum = 2;
		
		String output= "";//null;//"";
		
		ArrayList<Tolerance> tols = new ArrayList<Tolerance> ();
		Tolerance pmtol = new Tolerance(20.0f, true);//new Tolerance(20f, true);
		int MaxSpecNum = 1500;

		boolean accuracyFirst = false; 
		int minLength = 5, maxLength = 30, MaxPepLength = 20;
		int minNumPerLength = 0;
		
	
		ArrayList<String> paras = new ArrayList<String>();
		
		if(output != null){
			if(frag.length == 1)
				output = inputmgf[0] + ".out"+numPerLength+"_"+ charge;
			else{
				output = inputmgf[0] + "paired.out"+numPerLength+"_"+ charge;
			}
			if(setAccuracyThreshold>0){
				if(output != null) output = output + "_" + setAccuracyThreshold;
				maxGapNum = 200;
			}
		}
		
		//tols.add(new Tolerance(0.5f, false));
		
		for(String fr : frag){
			if(!fr.equals("HCD"))
				tols.add(new Tolerance(0.5f, false));
			else
				tols.add(new Tolerance(20f, true));
		}
		if(enz.equals("LysC")) enzyme = Enzyme.LysC;
		else if(enz.equals("AspN")) enzyme = Enzyme.AspN;
		else if(enz.equals("LysN")){
			enzyme = Enzyme.LysN;
			enz = "Trypsin";
		}
		
		if(numPerLength>20) 
			enzyme = null; 
	
		for(int i=0;i<frag.length;i++){
			String f = frag[i];
			String suffix = "_train.parsnewshort3";
			
			if(!f.equals("HCD")){
				if(enzyme.equals(Enzyme.LysC)){
					suffix = "LysC"+suffix;
				}else if(enzyme.equals(Enzyme.AspN)){
					suffix = "AspN"+suffix;
				}
				else{
					suffix = "Trypsin"+suffix;
				}
				//if(tols.get(i).getToleranceAsDa(500f)<0.02f) suffix = "FT"+suffix;
			}
			
			paras.add(System.getProperty("user.home") +  "/Dropbox/UniNovo/Pars" + System.getProperty("file.separator") + f+suffix);
		}
	
		//enzyme = null; // TODO
		PrintStream outfile = null;
		if(output != null) outfile = new PrintStream(output);
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		//aaSet = AminoAcidSet.getAminoAcidSetFromModFile("/home/kwj/workspace/MSGFDB/ModsPhospho.txt");
		/*aaSet = AminoAcidSet.getAminoAcidSetFromModFile("/home/kwj/workspace/MSGFDB/Mods.txt");
		//aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCysWithTerm();
		for(AminoAcid aa :  aaSet) System.out.println(aa + "\t" + aa.isModified() + "\t");
		for(ModifiedAminoAcid aa :  aaSet.getNTermVariableMods()) System.out.println(aa + "*\t" + aa.isResidueSpecific() + "\t");
		for(ModifiedAminoAcid aa :  aaSet.getNTermFixedMods()) System.out.println(aa + "**\t" + aa.isResidueSpecific() + "\t");
		
		for(Peptide pep : new Peptide("KLDELIEQ").getModifiedPeptides(aaSet))
			System.out.println(pep.getMass());
		
		System.out.println(new Peptide("KLDELIEQ").getMass());
		System.exit(1);
		*/
		
		float avgpeplen = 0;
		WindowFilter filter = new WindowFilter(6, 50);
		float correctCutLen = 0;
		float incorrectCutLen = 0;
		
		float[][] peakaccuracydivider = null;
		float[][] peakaccuracy = null;
		float[][] gpaccuracydivider = new float[50][10];
		float[][] gpcorrectaccuracy = new float[50][10];
		int[] numcorrect = new int[50];
		float[] avgsize = new float[50];
		float[] avglength = new float[50];
		float totalcorrectlength = 0;
		float[] avgIdellength = new float[50];
		float[] avgcorrectlength = new float[50];
		int[] sn = new int[50];
		int[] tsn = new int[50];
		int specnum = 0;
		float[] gapnum = new float[50];
		int qualifiedSpecnum = 0;
		int totalnumcorrect = 0;
		float[][] nodeaccuracydivider = new float[50][10];
		float[][] nodeaccuracy = new float[50][10];
		float[][] edgeaccuracydivider = new float[50][10];
		float[][] edgeaccuracy = new float[50][10];
		
		File[] specFiles = new File[frag.length];
		
		for(int i=0; i<frag.length; i++){
			specFiles[i] = new File(inputmgf[i]);
		}

		double time = 0;
		
		ArrayList<Iterator<Spectrum>> iterators = new ArrayList<Iterator<Spectrum>>();
		ArrayList<IPDGenerator> ipgs = new ArrayList<IPDGenerator>();
		
		for(int i=0; i< specFiles.length; i++){
			if(specFiles[i].isDirectory()) continue;
			System.out.println(specFiles[i]);
			ipgs.add(new IPDGenerator(paras.get(i), aaSet, i).filter(filter));
			Iterator<Spectrum> iterator = new SpectraIterator(specFiles[i].getAbsolutePath(), new MgfSpectrumParser());
			iterators.add(iterator);
		}
		
		
		if(charge > 0){
			peakaccuracydivider = new float[ipgs.get(0).getSigIonsOrderedByIntensityWithOutNoiseIon(charge).size()][10];
			peakaccuracy = new float[ipgs.get(0).getSigIonsOrderedByIntensityWithOutNoiseIon(charge).size()][10];
		}
		
		
		
		HashSet<String> peps = new HashSet<String>();
		if(accuracyFirst) SpectrumGraph.setAccuracyFirst();
		
		
		while(iterators.get(0).hasNext()){
			///////////////////////////////
			double etime = System.nanoTime(); 
			ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
			
			for(Iterator<Spectrum> iterator : iterators){
				spectra.add(iterator.next());
			
			}

			
			//System.out.println("****");
			if(charge > 0 && spectra.get(0).getCharge() != charge) continue;
			/*if(specnum < 500){
				if(spectra.get(0).getCharge() != 2) continue;
			}else
				if(spectra.get(0).getCharge() != 3) continue;
			*/
	
			Peptide annotation = spectra.get(0).getAnnotation();
		//	if(annotation == null) annotation = spectra.get(1).getAnnotation();
			if(annotation.isModified()) continue; //TODO recover
			
		//	annotation = new Peptide(spectra.get(0).getTitle());// TODO kill
		//	spectra.get(0).correctParentMass(annotation.getParentMass());
			
			System.out.println(annotation);
		//	if(spectra.size()>1) System.out.println(spectra.get(0).getAnnotationStr() + "\t" + spectra.get(1).getAnnotationStr());
			assert(spectra.size() == 1 || spectra.get(0).getParentMass() == spectra.get(1).getParentMass());
		
			
			SpectrumGraph graph = new SpectrumGraph(spectra, ipgs, enzyme, tols, pmtol, maxPTMNum > 0, true, true);
		//	System.out.println(graph.getNodes().size());
			UniNovo mg = new UniNovo(graph);
		
			ArrayList<DenovoReconstruction> can =  mg.getCandidateDenovoReconstructions(Math.max(numPerLength * 3, 100), minLength, maxLength, maxGapNum, maxPTMNum, accuracyFirst);
			ArrayList<DenovoReconstruction> out = new ArrayList<DenovoReconstruction>();
			
			float setaccuracy = mg.selectOutputFromCandidates(out, can, numPerLength, minNumPerLength, setAccuracyThreshold, maxPTMNum, accuracyFirst);
			time += System.nanoTime() - etime; 
						
			///////////////////////////////
		
			avgpeplen += annotation.size();
			peps.add(annotation.toString());
			int pepLength = 0;// spectra.get(0).getAnnotation().size();
			if(pepLength > MaxPepLength)  pepLength = MaxPepLength;
			
			if(specnum > MaxSpecNum){
				 break;
			}
			tsn[pepLength]++; specnum++;
			
			//if(specnum > 1000) break;
			
			HashMap<Peak, HashMap<IonType, Float>>  profile = graph.getProfile();
			PeakGenerator pg = new PeakGenerator(spectra.get(0));
			if(spectra.get(0).getAnnotation() == null)
				pg = new PeakGenerator(spectra.get(1));
			
			if(charge > 0 && spectra.size() == 1){
				
				//SpectrumParameter spar = new SpectrumParameter(spectra.get(0), ipgs.get(0).getMinCharge(),  ipgs.get(0).getMaxCharge());
				
				for(Peak p : profile.keySet()){
					HashMap<IonType, Float> prof = profile.get(p);
					PeakParameter ppar = new PeakParameter(p, spectra.get(0), 0);
				//	System.out.println(prof);
					for(int i=0; i< peakaccuracydivider.length; i++){
						IonType ion = ipgs.get(0).getSigIonsOrderedByIntensityWithOutNoiseIon(spectra.get(0).getCharge(), ppar.getBasePeakPartitionNum()).get(i);
						
						int bin = (int)(prof.get(ion) * peakaccuracydivider[i].length);
						peakaccuracydivider[i][Math.min(bin, peakaccuracydivider[i].length-1)] ++;
						//if(pg.isExplainedBy(p, ion, tols.get(0), pmtol)){
						//	peakaccuracy[i][Math.min(bin, peakaccuracy[i].length-1)] ++; //TODO recover
						//}
						/*else if(i==0 && prof.get(ion)>0.9f){
							for(IonType ion2 :  ipgs.get(0).getSigIonsOrderedByIntensityWithOutNoiseIon(charge)){
								if(pg.isExplainedBy(p, ion2, tols.get(0), pmtol)){
									System.out.println(ion2 + "****************");
								}
							}
						}*/
					}
				}
			}
			
			HashSet<Node> cNodes = graph.getCorrectNodes(annotation);
			HashSet<Edge> cEdges = graph.getCorrectEdges(annotation);
			
			for(int i=1; i<graph.getNodes().size()-1;i++){ // except sink source
				
				Node m = graph.getNodes().get(i);
				//if(!m.isEnzymatic()) continue;
				
				float accuracy = m.getAccuracy();
				int bin = (int)(accuracy * nodeaccuracydivider[pepLength].length);
				
				nodeaccuracydivider[pepLength][Math.min(bin, nodeaccuracydivider[pepLength].length-1)]++;

				int bin2 = 0 ;
				boolean nodeCorrect = cNodes.contains(m);
				
				for(Edge e : graph.getLinkedEdges(m)){
					//if(!e.getLeftNode().isCorrect) continue;
					
					bin2 =(int)(e.getAccuracy() * edgeaccuracydivider[pepLength].length);
					edgeaccuracydivider[pepLength][Math.min(bin2, edgeaccuracydivider[pepLength].length-1)]++;
					if(cEdges.contains(e)){
						edgeaccuracy[pepLength][Math.min(bin2, edgeaccuracydivider[pepLength].length-1)]++;
					}
				}
	
				if(nodeCorrect){	
					nodeaccuracy[pepLength][Math.min(bin, nodeaccuracydivider[pepLength].length-1)]++;
				}else{
			//		if(bin == nodeaccuracydivider.length) System.out.println("***" + accuracy + "\t" + cgp.toBitSetString() + "\t" + m + "\t" +(spec.getAnnotation().getNominalMass()-m.getNominalMass() ));
				}
				
			}
			
			int maxcorrectgplength = 0;

			ArrayList<Float> prms = new ArrayList<Float>();
			System.out.println(annotation);
			if(annotation.size() < 5){
				System.out.println("OOPS Short annotation!");
				System.exit(0);
			}
			float prm = 0;
			for(int l=0; l<annotation.size()-1; l++){
				prm += annotation.get(l).getMass();
				prms.add(prm);
			}
			float cutlen = prms.size();
			
		//	System.out.println(annotation.getModifiedPeptides(aaSet));
			for(DenovoReconstruction gp : out){
			//	System.out.println(gp.getGapMassRepresentation(aaSet));
				for(Node node : gp.getNodes()){
					if(node.getMass() == 0 || node.isSink()) continue;
					
					int i = Collections.binarySearch(prms, node.getMass());
					boolean c = false;
					
					if(i>=0){
						prms.remove(i);
						c = true;
					}
					else{
						i=-i-1;
						
						for(int j=-1;j<=0;j++){
							if(i+j>=0 && i+j < prms.size())
								if(Math.abs(prms.get(i+j) - node.getMass()) < node.getTol().getToleranceAsDa(Math.max(node.getMass(), gp.getPeptideMass() -node.getMass()))){
									prms.remove(i+j);
									c = true;
									break;
								}
						}
					}
					if(!c) incorrectCutLen++;
				}
				
				int bin = (int)(gp.getAccuracy()*gpaccuracydivider[0].length);
				
				gpaccuracydivider[pepLength][Math.min(bin, gpaccuracydivider[0].length-1)]++;
				if(gp.isCorrect(annotation)){ // 
					maxcorrectgplength = maxcorrectgplength > gp.length()? maxcorrectgplength : gp.length();
					gpcorrectaccuracy[pepLength][Math.min(bin, gpaccuracydivider[0].length-1)]++;
				}
			}
			
			 cutlen -= prms.size();
			 correctCutLen += cutlen;
			 
			/*
			if(out.isEmpty()){
				for(float prm : spec.getAnnotation().getPRMMasses(true, 0)){
					System.out.println(prm+"*");
				}
				System.out.println(spec.getPeptideMass() + "*");
				for(NewNode node : graph.getNodes()){
					System.out.println(node + "\t" + node.getLR() + "\t" + node.isCorrect);
				}
				
				for(NewNode node : graph.getNodes()){
					if(node.isCorrect || node.isSink()){
						System.out.println(node );
						for(NewEdge e : graph.getEdges(node)){
							System.out.println("\t" + e);
						}
					}
				}
				
			//	System.exit(0);
			}*/
		
			/*if(out.isEmpty()){
				System.out.println(graph.getNodes());
				System.out.println(spec.getAnnotation());
			
				 graph.getLR(spec.getAnnotation(), 0f,0f, charge);
					System.exit(0);
			}*/
			
			//System.out.println(graph.getNodes());
			
			if(setaccuracy < setAccuracyThreshold){
				System.out.println(spectra.get(0).getScanNum() + "\tfiltered out");
			/*	System.out.println(spec.getPeptideMass() + "\t" + spec.getAnnotationStr() + "\n" + graph.getNodes());
				for(float prm : spec.getAnnotation().getPRMMasses(true, 0)){
					System.out.print(prm+" ");
				}
				System.out.println(spec.getPeptideMass() + "*");*/
				continue;
			}

			qualifiedSpecnum++;
			boolean isCorrect = false;
			sn[pepLength]++;
			
			
			System.out.println(spectra.get(0).getScanNum() + "\t" + annotation + "\t" + graph.getScore(annotation));
			//for(float prm : spec.getAnnotation().getPRMMasses(true, 0)) System.out.print(prm + " ");
		//	System.out.println(graph.toFileString());
			
		//	System.out.println("LR of the correct one : " + graph.getLR(spec.getAnnotation(), 0f,0f, charge));
			
			int size = 0;
			
		//	ArrayList<Peptide> fragd = getFragmentedPeptidesFrom(spec, ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(charge), tol, maxRank);
			//System.out.println(fragd);
			
			float correctlen = 0;
			for(DenovoReconstruction gp : out){
				System.out.print(gp + " : " + Math.round(gp.getAccuracy()*100)+"% ");
				
			//	System.out.println(gp.getMatches(sa, aaSet, true));
				
				if(gp.isCorrect(annotation)){//
					isCorrect = true;
					System.out.println(" true*\t" + gp.length() + "\t" + maxcorrectgplength + "\t" + gp.getLR());
					if(correctlen == 0) correctlen = gp.length();
					//			lengthLRaccuracy[gp.length()][Math.max(0, Math.round(gp.getLR()))]++;
				}else{
					System.out.println(" false\t" + gp.length() +"\t"+  maxcorrectgplength + "\t" +  gp.getLR());						
				}
			//	if(Math.round(gp.getLR()) <0 || Math.round(gp.getLR()) >= lengthLRdivider[gp.length()].length);
			//	else lengthLRdivider[gp.length()][Math.round(gp.getLR())]++;
				
				size++;
				gapnum[pepLength] += gp.getGapNum();
			//	Node prev = Node.getNode(0);//
			//	for(Node n : gp.getNodes()){
			//		System.out.println("+" + n+"\t" + n.getLR() + "\t" + Edge.getEdge(prev, n).getLR());
			//		prev = n;
			//	}
			}
			
//totalcorrectlength
			
			if(isCorrect){
				avglength[pepLength] += correctlen;
				totalcorrectlength += correctlen;
				numcorrect[pepLength]++;
				totalnumcorrect++;
				if(output != null)  outfile.println(spectra.get(0).getCharge()+"\t"+spectra.get(0).getAnnotationStr() + "\t" + correctlen);
			}
			avgsize[pepLength] += size;
			
			avgIdellength[pepLength] += 0 + 1;
			avgcorrectlength[pepLength] += maxcorrectgplength;
				
			System.out.println(size + " " + isCorrect +" " + setaccuracy*100f);
		//	if(spec.getScanNum() == 176) System.exit(0);
		//	System.exit(0);
		}
		
		if(output != null)  outfile.println("###\t" + correctCutLen + "\t" + incorrectCutLen + "\t" + qualifiedSpecnum + "\t" + totalnumcorrect);
		if(output != null)  outfile.close();
		
		System.out.println("stat = [");
		for(int i=0;i<sn.length;i++){
			if(sn[i] > 0){
				System.out.println(i+" " +(float) numcorrect[i]/sn[i] + " " +(float) avglength[i]/numcorrect[i]);
			} 
		}
		System.out.println("];");
	
		System.out.print("%");	
		for(int i=0;i<sn.length;i++){
			if(sn[i] > 0){
				System.out.print(i+" ");
			}
		}
		System.out.println();
		System.out.println("accuracy=[");
		for(int j=0; j<gpcorrectaccuracy[0].length;j++){
		//	System.out.println("%" + j*gpaccuracydivider[0].length+"-"+ ((j+1)*gpaccuracydivider[0].length) + "% : ");
			for(int i=0;i<sn.length;i++){
				if(sn[i] > 0 && gpaccuracydivider[i][j]>50){
					System.out.print(gpcorrectaccuracy[i][j]/gpaccuracydivider[i][j]+ " ");
				}else System.out.print(Float.NaN + " ");
			}
			System.out.println();
		}
		System.out.println("];");
		
		System.out.println();
		System.out.println("nodeaccuracy=[");
		for(int j=0;j<nodeaccuracydivider[0].length;j++){
			//System.out.print(j + " ");
			for(int k=0; k<nodeaccuracydivider.length;k++){
				System.out.print(nodeaccuracy[k][j]/nodeaccuracydivider[k][j] + " ");
			}
			System.out.println();
		}
		System.out.println("];");
		
		System.out.println();
		System.out.println("edgeaccuracy=[");
		for(int j=0;j<edgeaccuracydivider[0].length;j++){
			//System.out.print(j + " ");
			for(int k=0; k<edgeaccuracydivider.length;k++){
				System.out.print(edgeaccuracy[k][j]/edgeaccuracydivider[k][j] + " ");
			}
			System.out.println();
		}
		
		System.out.println("];");
	
		
		System.out.println();
		if(charge>0){
			//System.out.println("% " + ipgs.get(0).getSigIonsOrderedByIntensityWithOutNoiseIon(charge));
			System.out.println("peakaccuracy=[");
			for(int j=0;j<peakaccuracydivider[0].length;j++){
				//System.out.print(j + " ");
				for(int k=0; k<peakaccuracydivider.length;k++){
					System.out.print(peakaccuracy[k][j]/peakaccuracydivider[k][j] + " ");
				}
				System.out.println();
			}
			System.out.println("];");
		}
		
		System.out.println("spec: " + specnum + "\tqspec : " + qualifiedSpecnum);
		System.out.println("accuracy: " + ((float)totalnumcorrect/qualifiedSpecnum));
		System.out.println("total accuracy: " + ((float)totalnumcorrect/specnum));
		System.out.println("length: " + totalcorrectlength/totalnumcorrect + "\t" + avgpeplen/specnum);
		System.out.println("correct cut number: " + correctCutLen + "\t" + incorrectCutLen);
	//	System.out.println("gapnum: " + totalgapnum / totalSize);
		System.out.println(time/(float)specnum/1e9f+" sec/spec ");
	}
}
