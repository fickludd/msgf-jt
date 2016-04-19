package uninovoOld;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.NominalMassFactory;
import msgf.ScoredSpectrum;
import msgf.ScoredSpectrumSum;
import msgf.Tolerance;
import msscorer.NewRankScorer;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Enzyme;
import msutil.Matter;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;

public class ScoredSpectrumForMSGF<T extends Matter>  implements ScoredSpectrum<T>{

	private SpectrumGraph graph;
	private HashMap<T, Integer> nodeScoreTable;
	
	public ScoredSpectrumForMSGF(SpectrumGraph graph){
		this.graph = graph;
		this.nodeScoreTable = new HashMap<T, Integer>();
	}
	

	
	
	@Override
	public int getEdgeScore(T curNode, T prevNode, float edgeMass) {
		
	/*	if(nodeTable.containsKey(curNode) && nodeTable.containsKey(prevNode)){
			Edge e = graph.getEdge(nodeTable.get(prevNode), nodeTable.get(curNode));
			
			if(e == null) return 0;
			return Math.round(e.getLR() - e.getRightNode().getLR());
		}
		*/
		return 0;
		/*
		ArrayList<Node> c = graph.getNodesWithinTolerance(curNode.getMass(), new Tolerance(0.5f, false));
		ArrayList<Node> p = graph.getNodesWithinTolerance(prevNode.getMass(), new Tolerance(0.5f, false));
		
		if(c.isEmpty() || p.isEmpty()) return 0;
		
		//float pm = graph.getSinkNode().getMass();
		
		float score = 0;
		
		boolean isConnected = false;
		
		for(Node cn : c){
			boolean isSubConnected = false;
			for(Node pn : p){
				if(graph.getEdge(pn, cn) == null){
					continue;
				}
				isSubConnected = true;
				break;
			}
			
			isConnected |= isSubConnected;
			if(!isSubConnected) score -= 0;
			
		}
		
	//	if(!isConnected) return -1000;
		//System.out.println(score);
		return Math.round(score);*/
	}

	@Override
	public int getNodeScore(T prm, T srm) {
		if(nodeScoreTable.containsKey(prm)) return nodeScoreTable.get(prm);
		
		float score = -200;
		
		Node cnode = null;
		
		for(Node node : graph.getNodesWithinTolerance(prm.getMass(), new Tolerance(0.5f, false))){ // tol = 0.5Da
			if(score < node.getLR()){
				score = node.getLR();
				cnode = node;
			}
		}
		
		float nopeakLR = 0;

		for(int i=0; i<graph.getTypeNum(); i++){
			nopeakLR += Node.getNoPeakLR(i, graph.getCorrectedCharge(i), prm.getMass(), graph.getSinkNode().getMass());
		}
		
		if(cnode == null) score = nopeakLR;
		
		
		score = Math.max(score, nopeakLR);
		nodeScoreTable.put(prm, Math.round(score*5));
		
		
		return Math.round(score*5);
	}
	
	//me : -9.963864	msgf : -12.011144	mix : -11.210074
	
	public static void main(String[] args) throws IOException{

		String prefix = "ETDTrypsin";
		String specFile = "/home/kwj/Dropbox/Test/" + prefix + "_test_ranked_-1.mgf_10.mgf";//%me : -12.774833	msgf : -11.868218	mix : -12.559194
		
		//specFile = "/home/kwj/Dropbox/Test/" + prefix + "_test.mgf";//%me : -12.774833	msgf : -11.868218	mix : -12.559194
		String para = "/home/kwj/Dropbox/PAR/"+ prefix + "_train.parsnew";
		int charge =2;//%%me : -13.8832245	msgf : -13.400685	mix : -14.120797
		//String param = "/home/kwj/Dropbox/CID_MSGF.param";
		Enzyme enzyme = Enzyme.TRYPSIN;//Enzyme.LysC ;
		
		Tolerance tol = new Tolerance(0.5f, false);
		if(prefix.contains("HCD")) tol = new Tolerance(20f, true);
		Tolerance pmtol = new Tolerance(20f, true);
		
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		NominalMassFactory factory  = new NominalMassFactory(aaSet, enzyme, 50);
		
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		ActivationMethod at = ActivationMethod.CID;
		if(prefix.contains("ETD")) at = ActivationMethod.ETD;
		if(prefix.contains("HCD")) at = ActivationMethod.HCD;
		
		NewRankScorer scorer = NewScorerFactory.get(at, enzyme);
		int sn = 0;//new NewRankScorer(param);//
		
		
		ArrayList<IPDGenerator> ipgs = new ArrayList<IPDGenerator>();
		ipgs.add(new IPDGenerator(para,aaSet, 0).filter(new WindowFilter(6, 50)));
		
		ArrayList<Tolerance> tols = new ArrayList<Tolerance>();
		tols.add(tol);
	//	if(tol.getToleranceAsDa(100) > 0.1f)
		for(IPDGenerator ipg : ipgs)
			ipg.setMaxIterationNum(1);		
		
		
		float a1 = 0,a2 = 0,a3 = 0;
		System.out.println("sp=[");
		while(iterator.hasNext()){
			 Spectrum spec = iterator.next();
			 if(spec.getCharge() != charge) continue;
			 if(sn > 100) break;
			 //if(spec.getAnnotation().size() < 13) continue;
			 // spec.correctParentMass();
			 
			 Annotation annotation = new Annotation("R."+ spec.getAnnotationStr() + ".A", aaSet);
			 
			 ScoredSpectrum<NominalMass> ss3 =  scorer.getScoredSpectrum(spec);
				
			 AminoAcidGraph g3 =  new AminoAcidGraph(factory, spec.getParentMass(), ss3);
	
			 GeneratingFunction<NominalMass> gf3 = new GeneratingFunction<NominalMass>(g3).enzyme(enzyme);
			 gf3.doNotBacktrack().doNotCalcNumber();
			 
			 gf3.computeGeneratingFunction();
			 float p3 = gf3.getSpectralProbability(annotation);
		 
			// if(p3 > 1e-9f) continue;
			 sn++;
					 
			 a3 += Math.log10(p3);
			 ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
			 
			 spectra.add(spec);
			 SpectrumGraph  graph = new SpectrumGraph(spectra, ipgs, enzyme, tols, pmtol, true, false, false);
			 ScoredSpectrum<NominalMass> ss =  new ScoredSpectrumForMSGF<NominalMass>(graph);
			
			 AminoAcidGraph g =  new AminoAcidGraph(factory, spec.getParentMass(), ss);
	
			 GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(g).enzyme(enzyme);
			 gf.doNotBacktrack().doNotCalcNumber();
			 gf.computeGeneratingFunction();
			 float p1 = gf.getSpectralProbability(annotation);
			
			 a1 += Math.log10(p1);
			 
			 ArrayList<ScoredSpectrum<NominalMass>> sl = new ArrayList<ScoredSpectrum<NominalMass>>();
			 sl.add(ss); sl.add(ss3);
			 
			 ScoredSpectrumSum<NominalMass> ss2 = new ScoredSpectrumSum<NominalMass>(sl);
			 AminoAcidGraph g2 =  new AminoAcidGraph(factory, spec.getParentMass(), ss2);
				
			 GeneratingFunction<NominalMass> gf2 = new GeneratingFunction<NominalMass>(g2).enzyme(enzyme);
			 gf2.doNotBacktrack().doNotCalcNumber();
			 gf2.computeGeneratingFunction();
			 
			 float p2 = gf2.getSpectralProbability(annotation);
			
			 a2 += Math.log10(p2);
			 //System.out.println(sn + " " + spec.getAnnotationStr() +  " >>");
			 System.out.println(Math.log10(p1) + "\t" + Math.log10(p3));			 
		}
		System.out.println("];\n%me : " + a1/sn+ "\tmsgf : " + a3/sn + "\tmix : " + a2/sn);
		
	}
	

}
