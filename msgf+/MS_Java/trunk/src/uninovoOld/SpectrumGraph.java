package uninovoOld;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.NominalMassFactory;
import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Constants;
import msutil.Enzyme;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.Spectrum;
import msutil.WindowFilter;
import uninovoOld.train.PeakGenerator;

public class SpectrumGraph {
	
	static private boolean accuracyFirst = false;
	
	private ArrayList<Node> nodes;
	private ArrayList<Node> allNodes;
	private int[] charges;
	
	private HashMap<Node, ArrayList<Edge>> adjList;
	private HashMap<ArrayList<Node>, Edge> edgeMap;
	private AminoAcidSet aaSet;
	
	private boolean ptmEdgeAllowed = false;
	private Enzyme enzyme = null;
	private HashMap<Peak, HashMap<IonType, Float>> profile = null;
	private int typeNum = 1;
	static private int num13C = 0;
	
	static public void setAccuracyFirst() {accuracyFirst = true;}
	
	public int getTypeNum(){ return typeNum; }
	
	private Node mergeNodes(ArrayList<Float> closeMasses, ArrayList<IonType> ions, HashMap<Float, float[]> massIonProbs, HashMap<Float, float[]> massNoiseProbs, float pm, int charge, Tolerance tol, Tolerance pmtol, Enzyme enzyme, int typeIndex){
		float center = 0;
		Node node = null;
		boolean isEnzymaticNode = false;
		boolean isSink = false;
		int sigIonSize = ions.size();
		float maxProbPrefix = 0;
		float maxprobSuffix = 0;
		float bestPrefix = -1, bestSuffix = -1;
		float[] ionProb = new float[sigIonSize];
		//float[] noiseProb = new float[sigIonSize];
		boolean isEmpty = true;
		
		ArrayList<Float> enzymaticNodes = getEnzymaticMasses(pm, enzyme);
		
		for(float m : closeMasses){			
			
			float[] p = massIonProbs.get(m);
			//float[] n = massNoiseProbs.get(m);
			
			isEnzymaticNode |= enzymaticNodes.contains(m);
			isSink |= m == pm;
			for(int i=0; i<ions.size();i++){
				float prob = p[i];
				if(!isEnzymaticNode && prob == 0) continue;
				
				if(i == 0){
					ionProb[i] = prob;
			//		noiseProb[i] = n[i];
				}else{
					ionProb[i] = Math.max(ionProb[i], prob);
					//noiseProb[i] = Math.min(noiseProb[i], n[i]);
				}
				if(ionProb[i]>0.0f) isEmpty = false;
				
				if(ions.get(i) instanceof IonType.SuffixIon){
					if(prob >= maxprobSuffix){
						maxprobSuffix = prob;
						bestSuffix = m;
					}
				}else{
					if(prob >= maxProbPrefix){
						maxProbPrefix = prob;
						bestPrefix = m;
					}
				}
			}
			
		}
		
		if(isEmpty && !isEnzymaticNode) return null;
		
		if(closeMasses.get(0) < pm/2){ // prefix is better
			center = bestPrefix;
			if(center < 0) center = bestSuffix;
		}else{
			center = bestSuffix;
			if(center < 0) center = bestPrefix;
		}
		
		if(isSink) center = pm;
		//if(isEnzymaticNode) System.out.println(center + "**");
		
		if(center >= 0 && center <= pm){
			if((SpectrumGraphComponent.useBinning(tol) && !isSink))// || (isSink&& useBinningSink))
				center = Math.round(center * Constants.INTEGER_MASS_SCALER) /   Constants.INTEGER_MASS_SCALER;
			
			node = new Node(center, charge, ionProb, tol, pm, isEnzymaticNode, typeIndex);
			if(isSink){
				node.setSink(pmtol);
			}
		}
		
		return node;
	}
	
	public HashMap<Peak, HashMap<IonType, Float>> getProfile() { return profile; } 
	
	static public void setNum13C(int n){ num13C = n;}
	
	// Spectrum spec is NOT updated
	private ArrayList<Node> generateNodes(Spectrum spec, IPDGenerator ipe, Enzyme enzyme, Tolerance tol, Tolerance pmtol, int typeIndex){
		
		spec.setRanksOfPeaks();
		Spectrum s = spec.getCloneWithoutPeakList();
	
		int maxRank = IPDGenerator.getMaxRank(s);
		
		for(Peak p : spec){ 
			if(p.getRank() < maxRank) 
				s.add(p.clone());
		}

		ArrayList<Node> tnodes = new ArrayList<Node>();
		profile = ipe.getIPD(s, 100);
		ArrayList<IonType> sigIons = ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(getCorrectedCharge(typeIndex));
	
		HashMap<Float, float[]> massIonProbs = new HashMap<Float, float[]>();
		//HashMap<Float, float[]> massNoiseProbs = new HashMap<Float, float[]>();
		
		ArrayList<Float> masses = new ArrayList<Float>(); 
		
		float pm = s.getPeptideMass();
		
		for(Peak peak : profile.keySet()){
			HashMap<IonType, Float> prof = profile.get(peak);
			for(int i=0; i<sigIons.size();i++){
				IonType ion = sigIons.get(i);
				float mass = PeakGenerator.getPrefixMass(peak, ion, s);
				if(mass >= pm || mass <= 0) continue;
				masses.add(mass);
				float[] probs = massIonProbs.get(mass);
				//float[] noiseProbs = massNoiseProbs.get(mass);
				
				if(probs == null)
					probs = new float[sigIons.size()];

			//	if(noiseProbs == null)
			//		noiseProbs = new float[sigIons.size()];

				
				probs[i] += prof.get(ion);
				massIonProbs.put(mass, probs);
				
			//	noiseProbs[i] += prof.get(IonType.NOISE);
			//	massNoiseProbs.put(mass, noiseProbs);
			}
		}

		for(float m :getEnzymaticMasses(pm, enzyme)){ // for enzymatic nodes
			if(masses.contains(m)) continue;
			masses.add(m);
			massIonProbs.put(m, new float[sigIons.size()]);
		//	massNoiseProbs.put(m, new float[sigIons.size()]);
		}
		
		float[] sourceSinkProbs = new float[sigIons.size()]; // for source sink
		for(int i=0; i<sigIons.size();i++){sourceSinkProbs[i] = 1;}
		masses.add(0f);
		massIonProbs.put(0f, sourceSinkProbs);
		//massNoiseProbs.put(0f, new float[sigIons.size()]);
		masses.add(pm);
		massIonProbs.put(pm, sourceSinkProbs);
		//massNoiseProbs.put(pm, new float[sigIons.size()]);
		
		Collections.sort(masses);
		
		ArrayList<Float> closeMasses = new ArrayList<Float>();
		for(int i=0; i<masses.size(); i++){
			float mass = masses.get(i);
			
			float delta = tol.getToleranceAsDa(Math.max(mass, pm-mass));
			if(closeMasses.isEmpty())
				closeMasses.add(mass);
			else if(!SpectrumGraphComponent.useBinning(tol) && Math.abs(closeMasses.get(0) - mass) <= 2*delta)
				closeMasses.add(mass);
			else if(NominalMass.toNominalMass(mass) == NominalMass.toNominalMass(closeMasses.get(0)))
				closeMasses.add(mass);
			
			
			else{
				Node tn = mergeNodes(closeMasses, sigIons, massIonProbs, null, pm,  getCorrectedCharge(typeIndex), tol, pmtol,  enzyme, typeIndex);
		
				if(tn != null) tnodes.add(tn);
				closeMasses.clear();
				closeMasses.add(mass);
				
			}
			
		}

		Node tn = mergeNodes(closeMasses, sigIons, massIonProbs, null, pm,  getCorrectedCharge(typeIndex),  tol, pmtol, enzyme, typeIndex);
		if(tn != null && tn.getAccuracy() > 0) tnodes.add(tn);
		
		Collections.sort(tnodes);
		
		return tnodes;
	}
	
	private ArrayList<Node> selectHighLRScoreNodes(ArrayList<Node> tnodes, int numHighLRNodes){
		
		ArrayList<Node> nodes = new ArrayList<Node>();
		/*
		int numPRM = 2;//Math.min(10, (int) (2/tol.getToleranceAsDa(spec.getPeptideMass())));

		float window = 100;
	   for(int peakIndex = 0; peakIndex < tnodes.size(); peakIndex++) {
		      int rank = 1;
		      
		      Node thisPeak = tnodes.get(peakIndex);
		      float thisMass = thisPeak.getMass();
		      float thisInten = thisPeak.getLR();
		      
		      // move left
		      int prevIndex = peakIndex-1;
		      while(prevIndex >= 0) {
		    	 Node prevPeak = tnodes.get(prevIndex);
		        if(thisMass - prevPeak.getMass() > window)    break;
		        if(prevPeak.getLR() > thisInten)            rank++;
		        prevIndex--;
		      }

		      // move right
		      int nextIndex = peakIndex+1;
		      while(nextIndex < tnodes.size()) {
		        Node nextPeak = tnodes.get(nextIndex);
		        if(nextPeak.getMass() - thisMass > window)    break;
		        if(nextPeak.getLR() > thisInten)            rank++;
		        nextIndex++;
		      }
		    
		     if(rank <= numPRM)  
		    	 nodes.add(thisPeak);
		      
		  }
			
		*/
		
		int maxNodeNum = Math.round(tnodes.get(tnodes.size()-1).getMass()/121.6f * 2.5f);
		if(accuracyFirst) maxNodeNum = 20;
	
		maxNodeNum = Math.min(maxNodeNum, 30);
		
		if(numHighLRNodes > 0)
			maxNodeNum = Math.max(maxNodeNum, numHighLRNodes);
		
		maxNodeNum = Math.min(maxNodeNum, tnodes.size());
	
		Collections.sort(tnodes, Collections.reverseOrder(SpectrumGraphComponent.LRComparator.get()));
		
		float lr = tnodes.get(maxNodeNum-1).getLR();
		
		for(Node node : tnodes){
			//if(!nodes.contains(node))
				nodes.add(node);
			if(node.getLR() < lr || nodes.size() > maxNodeNum+5) break;
		}
		
		/*ArrayList<Float> lrs = new ArrayList<Float>();
		
		for(Node node : tnodes){
			lrs.add(node.getLR());
		}
		
		Collections.sort(lrs, Collections.reverseOrder());

		for(Node node : tnodes){

			if(node.getLR() >= lrs.get(Math.min(lrs.size()-1, maxNodeNum)))
				if(!nodes.contains(node))nodes.add(node);
		}*/
		Collections.sort(nodes);
		Collections.sort(tnodes);
		return nodes;
	}
	
	static public ArrayList<Float> getEnzymaticMasses(float peptideMass, Enzyme enzyme){
		ArrayList<Float> ret = new ArrayList<Float>();
		if(enzyme == null) return ret;
	
		for(AminoAcid aa : enzyme.getResidues()){
			if(enzyme.isCTerm()){// sink side
				ret.add(peptideMass - aa.getMass());
			}else{
				ret.add(aa.getMass());
			}
		}
		return ret;
	}
	
	//Only For Training!!!
	public SpectrumGraph(Spectrum spec, IPDGenerator ipg, Enzyme enzyme, Tolerance tol, Tolerance pmtol, boolean ptmEdgeAllowed, boolean generateEdge){
		typeNum = 1;
	
		if(enzyme != null) SpectrumGraphComponent.setEnzyme(enzyme);
		
		this.enzyme = enzyme;
		
		charges = new int[typeNum];
		charges[0] = spec.getCharge();
		
		this.ptmEdgeAllowed = ptmEdgeAllowed;
		aaSet = ipg.getAASet();
		
		ArrayList<Node> tnodes = generateNodes(spec, ipg, enzyme, tol, pmtol, 0);

		allNodes = tnodes;
				
		nodes = selectHighLRScoreNodes(tnodes, 0);
		
		if(generateEdge) updateAdjList(spec.getPeptideMass(), enzyme, true);
	}
	

	public SpectrumGraph(ArrayList<Spectrum> spectra, ArrayList<IPDGenerator> ipgs, Enzyme enzyme, ArrayList<Tolerance> tols, Tolerance pmtol, boolean ptmEdgeAllowed, boolean generateEdge, boolean correctPM){
		ArrayList<Node> tnodes = new ArrayList<Node>();
		aaSet = ipgs.get(0).getAASet();
		typeNum = spectra.size();
		charges = new int[typeNum];
		this.ptmEdgeAllowed = ptmEdgeAllowed;
		
		if(enzyme != null) SpectrumGraphComponent.setEnzyme(enzyme);
		
		if(correctPM) correctParentMass(spectra,  ipgs, tols, pmtol);
			 
		if(pmtol.getToleranceAsDa(500f) > 0.5f)
			pmtol = new Tolerance(0.5f, false);
		
		for(int i=0; i<spectra.size(); i++){
			Spectrum s = spectra.get(i);
			charges[i] = s.getCharge();
			IPDGenerator ipg = ipgs.get(i);
			Tolerance tol = tols.get(i);
			
			ArrayList<Node> tn = generateNodes(s, ipg, enzyme, tol, pmtol, i);
			
			tnodes.addAll(tn);

		}
		
		if(spectra.size()>1) sumNodes(tnodes);
		
		//tnodes.get(tnodes.size()-1).setSink(pmtol);
		
		allNodes = tnodes;
		nodes = selectHighLRScoreNodes(tnodes, 0);
		this.enzyme = enzyme;
		if(generateEdge) updateAdjList(spectra.get(0).getPeptideMass(), enzyme, true);
		
	}
	
	SpectrumGraph(ArrayList<Node> nodes, Enzyme enzyme, Tolerance pmtol, AminoAcidSet aaSet, int[] charges, int typeNum, boolean generateEdge){
		this.aaSet = aaSet;
		//this.nodes = nodes;
		this.allNodes = nodes;
		this.enzyme = enzyme;
		this.charges = charges;
		this.typeNum = typeNum;
		
		if(enzyme != null) SpectrumGraphComponent.setEnzyme(enzyme);
		
		Node sink = allNodes.get(allNodes.size()-1);
		sink.setSink(pmtol);
		
		this.nodes = selectHighLRScoreNodes(allNodes, 0);
		
		if(generateEdge) updateAdjList(sink.getMass(), enzyme, false);	
	}
	
	private void correctParentMass(ArrayList<Spectrum> spectra,  ArrayList<IPDGenerator> ipgs, ArrayList<Tolerance> tols, Tolerance pmtol){
		
		float origPm = spectra.get(0).getParentMass();
		float pmTolAsDa = pmtol.getToleranceAsDa(origPm);
		float updateStep = 0.5f;
		int iterationNum = (int)(pmTolAsDa/updateStep);
		
		if(pmTolAsDa <= 0.5f && num13C == 0) return;
		
		
		for(int k=0; k<spectra.size(); k++){
			Spectrum s = spectra.get(k);
			if(s.getCharge() == 0) s.setCharge(ipgs.get(k).getMinCharge());
			
			int charge = Math.max(s.getCharge(), ipgs.get(k).getMinCharge());
			
			charge = Math.min(charge, ipgs.get(k).getMaxCharge());
			
			if(charge != s.getCharge()){

				float parentMass= s.getParentMass();
				s.setCharge(charge);
				s.correctParentMass(parentMass);

			}
		}
		
		
		NominalMassFactory factory  = new NominalMassFactory(aaSet, enzyme, 50);
		IPDGenerator.setForPMCorrection(true);
		
		float bestPm = 0;
		float maxScore = -200;
		//float specProb = 1e-9f;
		
		for(int i=0; i<=num13C; i++){
			float pm = origPm - (float)(i * Composition.ISOTOPE);//
			for(int j=-iterationNum;j<=iterationNum;j++){
		
				float cpm = pm + j * updateStep;
				
				ArrayList<Spectrum> newSpectra = new ArrayList<Spectrum>();
				
				for(Spectrum s : spectra){
					newSpectra.add(s.getDeepClone());
				}
				
				for(int k=0; k<newSpectra.size(); k++){
					Spectrum s = newSpectra.get(k);
					s.correctParentMass(cpm);
				//	ipgs.get(k).getIPD(s, 1)
				}

				SpectrumGraph  graph = new SpectrumGraph(newSpectra, ipgs, enzyme, tols, pmtol, true, false, false);
				
				graph.getSinkNode().setSink(pmtol);
				
				float score = 0;
				float t = (float) Math.pow(.5, 1+i+Math.abs(j));
				score += (float)Math.log(t/(1-t));
				
				/*Collections.sort(graph.getNodes(), Collections.reverseOrder(LRComparator.get()));
				
				
				for(int k=0; k<5;k++){
					Node n = graph.getNodes().get(k);
					if(n.isSink() || n.getMass() == 0) continue;
					score += n.getLR();
				}
				*/
			//	System.out.println(cpm + "\t" + accuracy);
				ScoredSpectrum<NominalMass> ss =  new ScoredSpectrumForMSGF<NominalMass>(graph);
				AminoAcidGraph g =  new AminoAcidGraph(factory, cpm, ss);
				GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(g).enzyme(enzyme);
				gf.doNotBacktrack().doNotCalcNumber().doNotCalcProb();
				gf.computeGeneratingFunction();
				//score += gf.getThresholdScore(specProb);
				score += gf.getMaxScore();
				
				if(maxScore < score){
					//System.out.println(cpm + "\t" + score);
					maxScore = score;
					bestPm = cpm;
				}
			}
		}
		
		for(int k=0; k<spectra.size(); k++){
			Spectrum s = spectra.get(k);
			
			s.correctParentMass(bestPm);
		}
		
		IPDGenerator.setForPMCorrection(false);
		//System.out.println("original pm: " + origPm + "\tcorrected pm: " + bestPm);
	}
	
	
	static public SpectrumGraph getSpectrumGraphFromPRMSpectrum(Spectrum s, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet, Enzyme enzyme){
		ArrayList<Node> nodes = new ArrayList<Node>();
		float pm = s.getPeptideMass();
		
		Node source = new Node(0, tol, pm, 0);
		source.setLR(20);
		nodes.add(source);
		for(Peak p : s){
			if(p.getMz() == 0) continue;
		
			
			Node node = new Node(p.getMz(), tol, pm, 0 );
			if(node.isWithinTolerance(s.getPeptideMass(), s.getPeptideMass(), pmtol)) continue;
			node.setLR(p.getIntensity());
		
			nodes.add(node);

		}
		Node sink = new Node(pm, tol, pm, 0);
		sink.setSink(pmtol);
		sink.setLR(20);
		nodes.add(sink);
		
		int[] charges = new int[1];
		charges[0] = s.getCharge();
		
		return new SpectrumGraph(nodes, enzyme, pmtol, aaSet, charges, 1, true);
	}
	
	
	private void sumNodes(ArrayList<Node> nodes){
		Collections.sort(nodes);
		
		Node prev = null;
		for(int i=0; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Node surv = node, remov = null;
			if(prev != null){
				Tolerance tol1 = prev.getTol();
				Tolerance tol2 = node.getTol();
				Tolerance tols, toll;
				
				if(tol1.getToleranceAsDa(node.getMass()) <= tol2.getToleranceAsDa(node.getMass())){
					tols = tol1; toll = tol2;
				}else{
					tols = tol2; toll = tol1;
				}
				
				if(Math.abs(node.getMass() - prev.getMass()) <= toll.getToleranceAsDa(node.getMass())){
					if(tol1.getToleranceAsDa(100) == tol2.getToleranceAsDa(100)){
						if(prev.getLR() > node.getLR()){
							surv = prev;
							remov = node;
						}else{
							surv = node;
							remov = prev;
						}
					}else if(tol1 == tols){
						surv = prev;
						remov = node;
					}else{
						surv = node;
						remov = prev;
					}
					surv.setAccuracy(Math.max(surv.getAccuracy(), remov.getAccuracy()));
				//	surv.setAccuracy(1-(1-surv.getAccuracy()) * (1-remov.getAccuracy()));
					//surv.calulateLRfromAccuracy();
					surv.setLR(surv.getLR() + remov.getLR());
				//	surv.calulateLRfromAccuracy();
					
					nodes.remove(nodes.indexOf(remov));
					i--;
					//System.out.println(remov + "\t" + surv);
				}
				
			}
			
			prev = surv;
		}
		Collections.sort(nodes);
	}
	
	private void updateAdjList(float pm, Enzyme enzyme, boolean updateNodePosteriorAccuracy){
		adjList = new HashMap<Node, ArrayList<Edge>>();
		edgeMap = new HashMap<ArrayList<Node>, Edge>();
	
		adjList.put(nodes.get(0), new ArrayList<Edge>());
		
		for(Node r : nodes){
			ArrayList<Edge> connectingEdges = new ArrayList<Edge>();
			ArrayList<Edge> aaEdges = new ArrayList<Edge>();
			
			for(Node l : nodes){
				if(l.getMass() >= r.getMass()) break;
				if(!adjList.containsKey(l)) continue;
				
				Edge connectingEdge = new Edge(l, r, pm, aaSet);
			
				if(updateNodePosteriorAccuracy){
					r.setMinAA(connectingEdge.getMinAANum());
					l.setMinAA(connectingEdge.getMinAANum());
				}
				
				ArrayList<Node> key = new ArrayList<Node>();
				key.add(l); key.add(r);
				
				if(enzyme != null){
					if(enzyme.isCTerm() && r.isSink() && !connectingEdge.isEnzymatic()){
						continue;
					}else if(enzyme.isNTerm() && l.getMass() == 0 && !connectingEdge.isEnzymatic()){
						continue;
					}
				}			
				edgeMap.put(key, connectingEdge);
	
				if(connectingEdge.isValid()){// &&  connectingEdge.getMass() < Edge.maxEdgeMass
					if(ptmEdgeAllowed){
						connectingEdges.add(connectingEdge);
					}
					else if(!connectingEdge.isPTM()){
						connectingEdges.add(connectingEdge);
						if(!connectingEdge.isGap())
							aaEdges.add(connectingEdge);
					}
				}
			}
			
			if(!aaEdges.isEmpty()){
				//r.setLR(r.getLR()*2);
				//adjList.put(r, aaEdges);
			}
			
			if(!connectingEdges.isEmpty()){
				/*Edge maxAccuraycGap = null;
				float maxAccuracy = -1;
				ArrayList<Integer> gapIndexes = new ArrayList<Integer>();
				
				for(int i=0; i<connectingEdges.size(); i++){
					Edge e = connectingEdges.get(i);
					if(e.isGap()){
						gapIndexes.add(i);
						if(maxAccuracy < e.getAccuracy()){
							maxAccuracy = e.getAccuracy();
							maxAccuraycGap = e;
						}
						connectingEdges.remove(i--);
					}
				}
				if(maxAccuraycGap != null) connectingEdges.add(maxAccuraycGap);*/
				adjList.put(r, connectingEdges);
			}
		}
		
		//setLRandAccuracy()
		

		if(updateNodePosteriorAccuracy){
			for(int i=1; i<nodes.size()-1;i++){
				nodes.get(i).setPosteriorLRnAccuracy();
			}
		}
		//System.out.println("++++" + edgeMap.size());
		for(ArrayList<Node> key : edgeMap.keySet()){
			edgeMap.get(key).setLRandAccuracy();
		}
	}
	
	public Enzyme getEnzyme(){return enzyme;}
	
	public int getCorrectedCharge(int type) { return charges[type]; }
	
	public boolean isEdgeGenerated(){
		return adjList != null;
	}
	
	public AminoAcidSet getAminoAcidSet(){
		return aaSet;
	}
	
	public Edge getEdge(Node l, Node r){
		if(edgeMap == null) return null;
		ArrayList<Node> key = new ArrayList<Node>();
		key.add(l); key.add(r);
		return edgeMap.get(key);
	}
	
	public ArrayList<Node> getNodes(){return nodes;}
	
	public ArrayList<Edge> getLinkedEdges(Node node){
		ArrayList<Edge> edges;
		edges = adjList.get(node);
		
		if(edges == null) return new ArrayList<Edge>();
		else	return edges;
	}
	
	public Node getSinkNode() {
		return nodes.get(nodes.size()-1);
	}
	
	public Node getSourceNode() {
		return nodes.get(0);
	}
	
	public float getScore(Peptide peptide){
		
		HashSet<Node> cNodes = getCorrectNodes(peptide);
		
		float score = this.getSourceNode().getLR();
		//System.out.println(score);
		float prm = 0;
		for(int i=0; i<peptide.size(); i++){
			boolean exist = false;
			prm += peptide.get(i).getMass();
			for(Node node : this.getNodesWithinTolerance(prm, new Tolerance(0.5f, false))){ // tol = 0.5Da
				if(cNodes.contains(node)){
					score += node.getLR();
					exist = true;
				//	break;
				}
			}
		//	System.out.println(score);
			if(!exist){
				for(int j=0; j<this.getTypeNum(); j++){
					score += Node.getNoPeakLR(j, this.getCorrectedCharge(j), prm, peptide.getMass());
				}
			}
			
		}
		//System.out.println(score+"*");
		return score;
	}
	
	public ArrayList<Node> getNodesWithinTolerance(float mass){
		return getNodesWithinTolerance(mass, null);
	}
	
	public ArrayList<Node> getNodesWithinTolerance(float mass, Tolerance tol){
		ArrayList<Node> ret = new ArrayList<Node>();
		Node m = new Node(mass, null, this.getSinkNode().getMass(), 0);
		int i = Collections.binarySearch(allNodes, m);
		if(i<0) i = -i-1;
		
		for(int j=-1;j<=1;j++){
			if(i+j >=0 && i+j<allNodes.size()){
				if(allNodes.get(i+j).isWithinTolerance(mass, getSinkNode().getMass(), tol))
				{
					ret.add(allNodes.get(i+j));
				}
			}
		}
		return ret;
	}
	
	
	public String toFileString(){
		String s = "";
		for(Node node : allNodes){
			s += node.toFileString();
		}
		return s;
	}
	
	
	static private void printUsageAndExit(){
		System.out.println("Usage : java -jar -Xmx1000m SpectrumGraph.jar" +
				"\n -i inputmgf files seperated by :" +
				"\n -o outputmgf file" +
				"\n -t ion tolerances seperated by : (each ending in ppm or Da)" +
				"\n -pt precursor ion tolerance (ending in ppm or Da)" +
				"\n -f fragmentation methods seperated by : (CID/ETD/HCD)" +
				"\n -e enzyme applied (0: No enzyme specificity, 1: Trypsin, 2: LysC 3: AspN)" +
				"\n -m allow modification edges - if set, nodes can be anywhere regardless of amino acid masses. (0: not allowed (default) 1: allowed )");
		System.out.println("Example : " +
				"\n java -jar -Xmx1000m SpectrumGraph.jar -i test.mgf -o output.mgf -t 0.1Da -pt 20ppm -f CID -e 1" +
				"\n java -jar -Xmx1000m SpectrumGraph.jar -i testCID.mgf:testHCD.mgf -o output2.mgf -t 0.1Da:20ppm -pt 20ppm -f CID:HCD");
		System.out.println("Applicable fragmentations : " +
				"\n CID : CID" +
				"\n ETD : ETD" +
				"\n HCD : HCD" );				
		System.out.println("Enzyme : " +
				"\n 0: no enzyme" +
				"\n 1: Trypsin" +
				"\n 2: LysC" +
				"\n 3: AspN");
		System.exit(0);
	}
	
	
	
	public ArrayList<Node> getAllNodes() {
		return allNodes;
	}
	
	public HashSet<Node> getCorrectNodes(Peptide pep){
		HashSet<Node> ret = new HashSet<Node>();
		
		for(Node node : getAllNodes()){
			if(node.isCorrect(pep))
				ret.add(node);
		}
		return ret;
	}
	
	public HashSet<Edge> getCorrectEdges(Peptide pep){
		HashSet<Edge> ret = new HashSet<Edge>();
		
		for(ArrayList<Node> k : edgeMap.keySet()){
			Edge e = edgeMap.get(k);
			if(e.isCorrect(pep)) ret.add(e);
		}
		
		return ret;
		
	}
	
	public float getParentMass(){		
		return getSinkNode().getMass() + (float)Composition.H2O;
	}
	
static public void main(String[] args) throws IOException{
		
		if(args.length < 1) printUsageAndExit();
		
		ArrayList<String> specFiles = new ArrayList<String>();//
		ArrayList<String> paras = new ArrayList<String>();//
		ArrayList<Tolerance> tols = new ArrayList<Tolerance>();//
		Tolerance pmtol = null;//
		String output = "";
		Enzyme enzyme = Enzyme.TRYPSIN;
		boolean allowPTM = false;
		
		ArrayList<String> frags = new ArrayList<String>();
		int mode = -1;
		for(String arg : args){
			if(arg.startsWith("-i")) mode = 0;
			else if(arg.startsWith("-o")) mode = 1;
			else if(arg.startsWith("-t")) mode = 2;
			else if(arg.startsWith("-pt")) mode = 3;
			else if(arg.startsWith("-f")) mode = 4;
			else if(arg.startsWith("-e")) mode = 5;
			else if(arg.startsWith("-m")) mode = 6;
			else{
				if(mode < 0) continue;
				if(mode == 0){
					for(String file : arg.split(":"))
						specFiles.add(file);
				}else if(mode == 1){
					output = arg;
				}else if(mode == 2){
					for(String tol : arg.split(":")){
						tols.add(Tolerance.parseToleranceStr(tol));
					}
				}else if(mode == 3){
					pmtol = Tolerance.parseToleranceStr(arg);
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
					}else{
						enzyme = Enzyme.AspN;
					}
				}else if(mode == 6){
					if(Integer.parseInt(arg) == 0) allowPTM = false;
				}
			}
		}
		
		for(int i=0;i<frags.size();i++){
			String frag = frags.get(i);
			String suffix = "_train.parsnew";
			
			if(!frag.equals("HCD")){
				if(enzyme.equals(Enzyme.LysC)){
					suffix = "LysC"+suffix;
				}else if(enzyme.equals(Enzyme.AspN)){
					suffix = "AspN"+suffix;
				}
				else{
					suffix = "Trypsin"+suffix;
				}
				if(tols.get(i).getToleranceAsDa(500f)<0.02f) suffix = "FT"+suffix;
			}
			
			
			paras.add("Pars" + System.getProperty("file.separator") + frag+suffix);
		}
		
		ArrayList<IPDGenerator> ipgs = new ArrayList<IPDGenerator>();
		
		WindowFilter filter = new WindowFilter(6, 50);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();

		for(int i=0; i<specFiles.size(); i++){
			ipgs.add(new IPDGenerator(paras.get(i), aaSet, i).filter(filter));
		}
		
		ArrayList<Iterator<Spectrum>> iterators = new ArrayList<Iterator<Spectrum>>();
		
		for(String file : specFiles){
			iterators.add(UniNovo.getSpectralIterator(file));
		}
		
		
		PrintStream outFile = new PrintStream(output);
				
		while(iterators.get(0).hasNext()){
			ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
			
			for(int i=0; i<iterators.size(); i++){
				Iterator<Spectrum> iterator = iterators.get(i);
				Spectrum spec = iterator.next();				
				spectra.add(spec);
			}
			
			//if(!isApplicable) continue;
			//System.out.println(allowPTM);
			SpectrumGraph graph = new SpectrumGraph(spectra, ipgs, enzyme, tols, pmtol, allowPTM, true, true);
			
			ArrayList<Node> nodes = graph.getNodes();
			if(nodes == null || nodes.isEmpty()) continue;
			
			System.out.println(spectra.get(0).getScanNum() + "\t" + spectra.get(0).getTitle() + "\t" + spectra.get(0).getPrecursorPeak().getMz());
			
			Spectrum outspec = spectra.get(0).getCloneWithoutPeakList();
			
			for(Node node :nodes){
				if(graph.adjList.containsKey(node))
					outspec.add(new Peak(node.getMass(), node.getLR(), 1));
			}
			outspec.outputMgf(outFile, true);
			
		}
		outFile.close();
		
	}

}
