package hybridsearch;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.NominalMassFactory;
import msgf.ScoredSpectrum;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Composition;
import msutil.Enzyme;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MS2SpectrumParser;
import parser.MSGFDBParser;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import parser.PSM;
import parser.PSMList;
import parser.SLGFParser;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class MSGFLib {
	
	//private int numCandidate = 10;
	private String msgfResultFile;
	private String slgfResultFile;
	private String outFile;
	private String specDir;
	private SuffixArray sa;
	
	static private AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
	static private NominalMassFactory factory = new NominalMassFactory(aaSet, Enzyme.TRYPSIN, 200); //TODO opt later
	static private NewRankScorer scorer = NewScorerFactory.get(ActivationMethod.CID, Enzyme.TRYPSIN);
	
	static private boolean isRead = false;
	
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
	
	static protected void setAlpha(float a){
		HybridPvalue.setAlpha(a);
	}
	
	public MSGFLib(String outFile, String msgfResultFile, String fasta, String slgfResultFile, String paraFile, String specDir){
		this.outFile = outFile;
		this.msgfResultFile = msgfResultFile;
		this.slgfResultFile = slgfResultFile;
		this.specDir = specDir;
		
		this.sa = new SuffixArray(new SuffixArraySequence(fasta));
		
		read(paraFile);
	}
	
	public void run(){
		PSMList<PSM> msgfPSMs = MSGFDBParser.parse(msgfResultFile);
		PSMList<PSM> slgfPSMs = SLGFParser.parse(slgfResultFile);
		//PSMList<PSM> hybridPSMs = new PSMList<PSM>();
		
		HashMap<String, PSMList<PSM>> msgfPSMsPerSpec = new HashMap<String, PSMList<PSM>>();
		HashMap<String, PSMList<PSM>> slgfPSMsPerSpec = new HashMap<String, PSMList<PSM>>();
		HashSet<String> allSpecFiles = new HashSet<String>();
		HashMap<String, Spectrum> allSpecs = new HashMap<String, Spectrum>();
		
		StringBuffer b = new StringBuffer();
		b.append(specDir);
		b.append(System.getProperty("file.separator"));
		
		
		for(PSM psm : msgfPSMs){
			StringBuffer k = new StringBuffer(b);
			k.append(psm.getSpecFileName());
			allSpecFiles.add(k.toString());
			k.append(":"); 
			k.append(psm.getScanNum());
			
			String e = k.toString();
			
			if(!msgfPSMsPerSpec.containsKey(e))
				msgfPSMsPerSpec.put(e, new PSMList<PSM>());
			/*else{
				System.out.println(msgfPSMsPerSpec.get(e).get(0).getPrecursorMz());
				System.out.println(e);
				System.out.println(psm.getPrecursorMz());
				System.exit(1);
			}*/
			msgfPSMsPerSpec.get(e).add(psm);
			
		}
	
		for(PSM psm : slgfPSMs){
			StringBuffer k = new StringBuffer(b);
			k.append(psm.getSpecFileName()); 
			k.append(":"); 
			k.append(psm.getScanNum());

			String e = k.toString();
			if(!slgfPSMsPerSpec.containsKey(e))
				slgfPSMsPerSpec.put(e, new PSMList<PSM>());
			
			slgfPSMsPerSpec.get(e).add(psm);
		}
			
		for(String specFile : allSpecFiles){
			Iterator<Spectrum> iterator = getSpectralIterator(specFile);
			while(iterator.hasNext()){
				Spectrum s = iterator.next();
				StringBuffer e = new StringBuffer(specFile);
				e.append(":");
				e.append(s.getScanNum()); 
				allSpecs.put(e.toString(), s);
			}
		}
		
		System.out.println("Hashing Done : " + allSpecs.size());
		
		
		try {
			PrintStream ps = new PrintStream(outFile);
			ps.println("#SpecFileName\tScan#\tCharge\tPrecursorMz\tPeptide\tProtein\tMSGFpValue\tSLGFpValue\tpValue\tIsEPvalue");
			int i=0;
			for(String k : allSpecs.keySet()){
			//	System.out.println("\t" + i++);
				//if(i<6122) continue; 
			//	if(allSpecs.get(k).getScanNum() != 1953) continue;
					
				PSMList<PSM> hpsm = getPSMsWithPValues(msgfPSMsPerSpec.get(k), slgfPSMsPerSpec.get(k), allSpecs.get(k));
			//	System.out.println(hpsm.size());
				
				for(PSM psm : hpsm){
					ps.println(psm.getSpecFileName() + "\t" + psm.getScanNum() + "\t" + psm.getCharge() + "\t" + psm.getPrecursorMz()
							+"\t" + psm.getPeptide().toStringWithModification(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys()) + "\t" + psm.getProtein() + "\t" +  psm.getScore("M") + "\t" + psm.getScore("S") + "\t" + psm.getProbScore() + "\t" + psm.getScore("IsEPvalue"));
				}
			}	

			ps.close();
	
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	static protected void read(String file){ // read trained parameters
		if(isRead) return;
		//System.out.println("Parameter file read...");
		try {
			BufferedLineReader in = new BufferedLineReader(file);
			String s;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#")){
					setAlpha(Float.parseFloat(s.split("\t")[1]));
					continue;
				}
				ConditionalProbability cp = new ConditionalProbability(s);
				cp.set(cp.get());
			}
			
			in.close();
			isRead = true;
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private PSMList<PSM> getPSMsWithPValues(PSMList<PSM> msgfPSMs, PSMList<PSM> slgfPSMs, Spectrum spec){
		PSMList<PSM> psms = new PSMList<PSM>();
		HashMap<Peptide, PSM> msgfPeps = new HashMap<Peptide, PSM>();
		HashMap<Peptide, PSM> slgfPeps = new HashMap<Peptide, PSM>();
		HashSet<Peptide> allPeps = new HashSet<Peptide>();
		
		if(msgfPSMs != null){
			for(PSM psm : msgfPSMs){
				msgfPeps.put(psm.getPeptide(), psm);
			}
		}
		if(slgfPSMs != null){
			for(PSM psm : slgfPSMs){
				slgfPeps.put(psm.getPeptide(), psm);
			}
		}
		
		allPeps.addAll(msgfPeps.keySet());
		allPeps.addAll(slgfPeps.keySet());
		
	//	System.out.println(spec.getScanNum() + "\t" + msgfPSMs);
		float pm = spec.getParentMass(); 
		
		ScoredSpectrum<NominalMass> scoreSpec =  scorer.getScoredSpectrum(spec);
		AminoAcidGraph graph1 =  new AminoAcidGraph(factory, spec.getParentMass(), scoreSpec); // TODO parameterize
		spec.correctParentMass(pm - (float)Composition.ISOTOPE);
		scoreSpec =  scorer.getScoredSpectrum(spec);
		AminoAcidGraph graph2 =  new AminoAcidGraph(factory, spec.getParentMass(), scoreSpec);
		spec.correctParentMass(pm - (float)Composition.ISOTOPE2);
		scoreSpec =  scorer.getScoredSpectrum(spec);
		AminoAcidGraph graph3 =  new AminoAcidGraph(factory, spec.getParentMass(), scoreSpec);
		spec.correctParentMass(pm + (float)Composition.ISOTOPE);
		scoreSpec =  scorer.getScoredSpectrum(spec);
		AminoAcidGraph graph4 =  new AminoAcidGraph(factory, spec.getParentMass(), scoreSpec);
		spec.correctParentMass(pm + (float)Composition.ISOTOPE2);
		scoreSpec =  scorer.getScoredSpectrum(spec);
		AminoAcidGraph graph5 =  new AminoAcidGraph(factory, spec.getParentMass(), scoreSpec);

		GeneratingFunction<NominalMass> gf1 = new GeneratingFunction<NominalMass>(graph1).enzyme(Enzyme.TRYPSIN);
		gf1.doNotBacktrack().doNotCalcNumber().computeGeneratingFunction();
		GeneratingFunction<NominalMass> gf2 = new GeneratingFunction<NominalMass>(graph2).enzyme(Enzyme.TRYPSIN);
		gf2.doNotBacktrack().doNotCalcNumber().computeGeneratingFunction();
		GeneratingFunction<NominalMass> gf3 = new GeneratingFunction<NominalMass>(graph3).enzyme(Enzyme.TRYPSIN);
		gf3.doNotBacktrack().doNotCalcNumber().computeGeneratingFunction();
		GeneratingFunction<NominalMass> gf4 = new GeneratingFunction<NominalMass>(graph4).enzyme(Enzyme.TRYPSIN);
		gf4.doNotBacktrack().doNotCalcNumber().computeGeneratingFunction();
		GeneratingFunction<NominalMass> gf5 = new GeneratingFunction<NominalMass>(graph5).enzyme(Enzyme.TRYPSIN);
		gf5.doNotBacktrack().doNotCalcNumber().computeGeneratingFunction();
		
		for(Peptide pep : allPeps){
			PSM slgfpsm = slgfPeps.get(pep);
			PSM msgfpsm = msgfPeps.get(pep);
			PSM psm = null;
			HybridPvalue pval = null;
			boolean isExpected = false;
		
			
			if(slgfpsm == null){
			//	System.out.println(pep + "\t" + msgfpsm.getProbScore());
				pval = new HybridPvalue(msgfpsm.getProbScore());
				psm = msgfpsm;
				isExpected = true;
				psm.score("M", msgfpsm.getProbScore());
				psm.score("S", -1);
			//	System.out.print("D");
			}else if(msgfpsm == null){
				//System.out.print(2);

				Annotation annotation = new Annotation(null, slgfpsm.getPeptide() , null);
		
				float specProb = gf1.getSpectralProbability(annotation);
				specProb = Math.min(specProb, gf2.getSpectralProbability(annotation));
				specProb = Math.min(specProb, gf3.getSpectralProbability(annotation));
				specProb = Math.min(specProb, gf4.getSpectralProbability(annotation));
				specProb = Math.min(specProb, gf5.getSpectralProbability(annotation));
				
				pval = new HybridPvalue(specProb, slgfpsm.getProbScore());
				psm = slgfpsm;
			//	ArrayList<String> matchedProtSeq = sa.getAllMatchingAnnotations(slgfpsm.getPeptideStr().toUpperCase());
				if(slgfpsm.getScore("T/D") == 1)
					psm.protein("REV");
				
				psm.score("M", specProb);
				psm.score("S", slgfpsm.getProbScore());
				//	System.out.println(pep.getParentMass() + "\t" + spec.getParentMass() + "\t" +annotation + "\t"+ specProb + "\t" + slgfpsm.getProbScore());
				//NQAALNPR
				//else{
				//	for(String pro : matchedProtSeq){
				//		psm.protein(pro);
				//		break;
				//	}
				//}
			//	System.out.print("D");
			}else{
				//System.out.print(3);
	//			System.out.println(pep + "\t" + msgfpsm.getProbScore() + "\t" + slgfpsm.getProbScore());
				pval = new HybridPvalue(msgfpsm.getProbScore(), slgfpsm.getProbScore());
				psm = msgfpsm;
				psm.score("M", msgfpsm.getProbScore());
				psm.score("S", slgfpsm.getProbScore());
				//System.out.print("D");
			}
			
			float pv = pval.get();
			//System.out.println(isExpected + "\t" + pv + "\t" + psm.getProbScore());
			psm.probScore(pv);
			if(isExpected) psm.score("IsEPvalue", 1);
			else psm.score("IsEPvalue", 0);
			psms.add(psm);
			//System.out.println("Done");
		}
		
	//	System.out.println("Done");
		return psms;
	}
	
	
}
