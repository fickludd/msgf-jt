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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import uninovo.SpectrumGraphComponent.LRComparator;
import uninovo.util.AminoAcid;
import uninovo.util.AminoAcidSet;
import uninovo.util.Composition;
import uninovo.util.Constants;
import uninovo.util.Enzyme;
import uninovo.util.IonType;
import uninovo.util.NominalMass;
import uninovo.util.Peak;
import uninovo.util.Peptide;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;

/**
 * SpectrumGraph class represents a spectrum graph on which de novo sequences are generated
 * @author kyowon
 */

public class SpectrumGraph {

	/** The num13 c. */
	static private int num13C = 0;
	
	/**
	 * Filter out low scoring nodes
	 * @param tnodes nodes
	 * @return filtered nodes
	 */
	static private ArrayList<Node> filterNodes(ArrayList<Node> tnodes){
		
		ArrayList<Node> nodes = new ArrayList<Node>();
		int maxNodeNum = Math.round(tnodes.get(tnodes.size()-1).getMass()/121.6f * 2.5f);
		maxNodeNum = Math.min(maxNodeNum, 30);
	
		maxNodeNum = Math.min(maxNodeNum, tnodes.size());
		Collections.sort(tnodes, Collections.reverseOrder(SpectrumGraphComponent.LRComparator.get()));
		float lr = tnodes.get(maxNodeNum-1).getLR();
		for(Node node : tnodes){
				nodes.add(node);
			if(node.getLR() < lr || nodes.size() > maxNodeNum+5) break;
		}
		Collections.sort(nodes);
		Collections.sort(tnodes);
		return nodes;
	}
	
	/**
	 * Get prefix masses of enzymatic cleavage sites
	 * @param peptideMass mass of peptide
	 * @param enzyme enzyme
	 * @return masses 
	 */
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
	
	/**
		 * Get a node from close masses  
		 * @param closeMasses masses gathered within a given tolerance
		 * @param ions considered ion types
		 * @param fpvs FPVs for masses
		 * @param pm parent mass
		 * @param charge charge of spectrum
		 * @param tol MS2 tolerance
		 * @param pmtol MS1 tolerance
		 * @param enzyme enzyme
		 * @param fragmentationMethodIndex fragmentation method index
		 * @return a new node 
		 */
	static private Node getOneNode(ArrayList<Float> closeMasses, ArrayList<IonType> ions, HashMap<Float, float[]> fpvs, float pm, int charge, Tolerance tol, Tolerance pmtol, Enzyme enzyme, int fragmentationMethodIndex){
		float center = 0;
		Node node = null;
		boolean isEnzymaticNode = false;
		boolean isSink = false;
		int sigIonSize = ions.size();
		float maxProbPrefix = 0;
		float maxprobSuffix = 0;
		float bestPrefix = -1, bestSuffix = -1;
		float[] fpv = new float[sigIonSize];
		boolean isEmpty = true;
		
		ArrayList<Float> enzymaticNodes = getEnzymaticMasses(pm, enzyme);		
		for(float m : closeMasses){						
			float[] p = fpvs.get(m);			
			isEnzymaticNode |= enzymaticNodes.contains(m);
			isSink |= m == pm;
			for(int i=0; i<ions.size();i++){
				float prob = p[i];
				if(!isEnzymaticNode && prob == 0) continue;
				
				if(i == 0){
					fpv[i] = prob;
				}else{
					fpv[i] = Math.max(fpv[i], prob);
				}
				if(fpv[i]>0.0f) isEmpty = false;
				
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
		if(center >= 0 && center <= pm){
			if((SpectrumGraphComponent.useBinning(tol) && !isSink))
				center = Math.round(center * Constants.INTEGER_MASS_SCALER) /   Constants.INTEGER_MASS_SCALER;
			
			node = new Node(center, charge, fpv, tol, pm, isEnzymaticNode, fragmentationMethodIndex);
			if(isSink){
				node.setSink(pmtol);
			}
		}
		return node;
	}
	
	/**
	 * Set the number of considered 13C 
	 * @param n number
	 */
	static public void setNum13C(int n){ num13C = n;}
	
	/** The aa set. */
	private AminoAcidSet aaSet;
	
	/** The adj list. */
	private HashMap<Node, ArrayList<Edge>> adjList; 
	
	/** The all nodes. */
	private ArrayList<Node> allNodes;

	/** The charges. */
	private int[] charges;
	
	/** The edge map. */
	private HashMap<ArrayList<Node>, Edge> edgeMap;
	   
	   /** The filtered nodes. */
	private ArrayList<Node> filteredNodes;	
	
	// edges representing PTMs - for future use
	/** The ptm edge allowed. */
	private boolean ptmEdgeAllowed = false;
	
	/**
	 * Constructor 
	 * @param spectra spectra possibly from different fragmentation methods
	 * @param pipgs peak ion prob generators to score nodes on the spectrum graph
	 * @param enzyme enzyme
	 * @param tols MS2 tolerances
	 * @param pmtol MS1 tolerance
	 * @param ptmEdgeAllowed edges representing PTMs are allowed? - for future use
	 * @param generateEdge decides if edges will be generated or not
	 * @param correctPM decides if parent mass should be corrected
	 */
	 public SpectrumGraph(ArrayList<Spectrum> spectra, ArrayList<PeakIonProbabilityGenerator> pipgs, Enzyme enzyme, ArrayList<Tolerance> tols, Tolerance pmtol, boolean ptmEdgeAllowed, boolean generateEdge, boolean correctPM){
		ArrayList<Node> tnodes = new ArrayList<Node>();
		aaSet = pipgs.get(0).getAASet();
		charges = new int[spectra.size()];
		this.ptmEdgeAllowed = ptmEdgeAllowed;		
		if(enzyme != null) SpectrumGraphComponent.setEnzyme(enzyme);		
		if(correctPM) correctParentMass(spectra,  pipgs, enzyme, tols, pmtol);			 
		if(pmtol.getToleranceAsDa(500f) > 0.5f)
			pmtol = new Tolerance(0.5f, false);		
		for(int i=0; i<spectra.size(); i++){
			Spectrum s = spectra.get(i);
			charges[i] = s.getCharge();
			PeakIonProbabilityGenerator pipg = pipgs.get(i);
			Tolerance tol = tols.get(i);			
			ArrayList<Node> tn = generateNodes(s, pipg, enzyme, tol, pmtol, i);
			tnodes.addAll(tn);
		}
		
		if(spectra.size()>1) mergeNodes(tnodes);		
		allNodes = tnodes;
		
		filteredNodes = filterNodes(tnodes);
		if(generateEdge) updateAdjList(enzyme, true);
	}
	
	/**
	 * Constructor for training
	 * @param spec spectrum
	 * @param pipg peak ion prob generator to score nodes on the spectrum graph
	 * @param enzyme enzyme
	 * @param tol MS2 tolerance
	 * @param pmtol MS1 tolerance
	 * @param ptmEdgeAllowed edges representing PTMs are allowed? - for future use
	 * @param generateEdge decides if edges will be generated or not
	 */
	SpectrumGraph(Spectrum spec, PeakIonProbabilityGenerator pipg, Enzyme enzyme, Tolerance tol, Tolerance pmtol, boolean ptmEdgeAllowed, boolean generateEdge){	
		if(enzyme != null) SpectrumGraphComponent.setEnzyme(enzyme);			
		charges = new int[1];
		charges[0] = spec.getCharge();		
		this.ptmEdgeAllowed = ptmEdgeAllowed;
		aaSet = pipg.getAASet();		
		ArrayList<Node> tnodes = generateNodes(spec, pipg, enzyme, tol, pmtol, 0);
		allNodes = tnodes;				
		filteredNodes = filterNodes(tnodes);		
		if(generateEdge) updateAdjList(enzyme, true);
	}
	
	/**
	 * Correct parent mass 
	 * @param spectra spectra 
	 * @param pipgs peak ion prob generators
	 * @param enzyme enzyme
	 * @param tols MS2 tolerances
	 * @param pmtol MS1 tolerance 
	 */
	private void correctParentMass(ArrayList<Spectrum> spectra,  ArrayList<PeakIonProbabilityGenerator> pipgs, Enzyme enzyme, ArrayList<Tolerance> tols, Tolerance pmtol){
		
		float origPm = spectra.get(0).getParentMass();
		float pmTolAsDa = pmtol.getToleranceAsDa(origPm);
		float updateStep = 0.5f;
		int iterationNum = (int)(pmTolAsDa/updateStep);
		
		if(pmTolAsDa <= 0.5f && num13C == 0) return;
		
		
		for(int k=0; k<spectra.size(); k++){
			Spectrum s = spectra.get(k);
			if(s.getCharge() == 0) s.setCharge(pipgs.get(k).getMinCharge());	
			int charge = Math.max(s.getCharge(), pipgs.get(k).getMinCharge());
			charge = Math.min(charge, pipgs.get(k).getMaxCharge());
			
			if(charge != s.getCharge()){
				float parentMass= s.getParentMass();
				s.setCharge(charge);
				s.correctParentMass(parentMass);
			}
		}
		
		
		PeakIonProbabilityGenerator.setForPMCorrection(true);
		
		float bestPm = origPm;
		float maxScore = -200;
		//float isotopeProb = 0.9f;
		//float zeroUpdateProb = 0.9f;
		
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
				}

				SpectrumGraph  graph = new SpectrumGraph(newSpectra, pipgs, enzyme, tols, pmtol, false, true, false);
				graph.getSinkNode().setSink(pmtol);
				
				float score = 0;
				
				float t = (float) Math.pow(.5, 1+i+Math.abs(j));
				score += (float)Math.log(t/(1-t));
				
				Collections.sort(graph.getFilteredNodes(), Collections.reverseOrder(LRComparator.get()));
				for(int k=0; k<5;k++){
					Node n = graph.getFilteredNodes().get(k);
					if(n.isSink() || n.getMass() == 0) continue;
					score += n.getLR();
				}
		
				if(maxScore < score){
					maxScore = score;
					bestPm = cpm;
				}
			}
		}
		
		for(int k=0; k<spectra.size(); k++){
			Spectrum s = spectra.get(k);			
			s.correctParentMass(bestPm);
		}
		System.out.println("\t PM corrected from " + origPm + " to " + bestPm);
		PeakIonProbabilityGenerator.setForPMCorrection(false);
	}
	
	/**
	 * Generate nodes of this spectrum graph.
	 *
	 * @param spec spectrum
	 * @param pipg peak ion prob generator to score nodes on the spectrum graph
	 * @param enzyme enzyme
	 * @param tol the tol
	 * @param pmtol MS1 tolerance
	 * @param fragmentationMethodIndex index for different fragmentation methods
	 * @return the node set
	 */
	private ArrayList<Node> generateNodes(Spectrum spec, PeakIonProbabilityGenerator pipg, Enzyme enzyme, Tolerance tol, Tolerance pmtol, int fragmentationMethodIndex){
		
		spec.setRanksOfPeaks();
		Spectrum s = spec.getCloneWithoutPeakList();
		int maxRank = PeakIonProbabilityGenerator.getMaxRank(s);
		for(Peak p : spec){ 
			if(p.getRank() < maxRank) 
				s.add(p.clone());
		}

		ArrayList<Node> tnodes = new ArrayList<Node>();
		HashMap<Peak, HashMap<IonType, Float>> peakIonProbMap = pipg.getPIPs(s);
		ArrayList<IonType> sigIons = pipg.getSigIonsOrderedByIntensityWithOutNoiseIon(getBoundedCharge(fragmentationMethodIndex));
	
		HashMap<Float, float[]> fpvs = new HashMap<Float, float[]>();		
		ArrayList<Float> masses = new ArrayList<Float>(); 
		
		float pm = s.getPeptideMass();
		
		for(Peak peak : peakIonProbMap.keySet()){
			HashMap<IonType, Float> ionProbs = peakIonProbMap.get(peak);
			for(int i=0; i<sigIons.size();i++){
				IonType ion = sigIons.get(i);
				float mass = PeakGenerator.getPrefixMass(peak, ion, s);
				if(mass >= pm || mass <= 0) continue;
				masses.add(mass);
				float[] fpv = fpvs.get(mass);
				if(fpv == null)
					fpv = new float[sigIons.size()];
				fpv[i] += ionProbs.get(ion);
				fpvs.put(mass, fpv);
			
			}
		}
		for(float m :getEnzymaticMasses(pm, enzyme)){ // for enzymatic nodes
			if(masses.contains(m)) continue;
			masses.add(m);
			fpvs.put(m, new float[sigIons.size()]);
		}
		
		float[] sourceSinkProbs = new float[sigIons.size()]; // for source sink
		for(int i=0; i<sigIons.size();i++){sourceSinkProbs[i] = 1;}
		masses.add(0f);
		fpvs.put(0f, sourceSinkProbs);		
		masses.add(pm);
		fpvs.put(pm, sourceSinkProbs);
		
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
				Node tn = getOneNode(closeMasses, sigIons, fpvs, pm,  getBoundedCharge(fragmentationMethodIndex), tol, pmtol,  enzyme, fragmentationMethodIndex);
				if(tn != null) tnodes.add(tn);
				closeMasses.clear();
				closeMasses.add(mass);				
			}			
		}

		Node tn = getOneNode(closeMasses, sigIons, fpvs, pm,  getBoundedCharge(fragmentationMethodIndex),  tol, pmtol, enzyme, fragmentationMethodIndex);
		if(tn != null && tn.getAccuracy() > 0) tnodes.add(tn);
		Collections.sort(tnodes);		
		return tnodes;
	}
	
	/**
	 * Return all nodes (unfiltered)
	 * @return all nodes
	 */
	public ArrayList<Node> getAllNodes() {
		return allNodes;
	}
	
	/**
	 * Get bounded charge values - bounds are specified by trained parameters; for example, if only charge 2 
	 * spectra are learned, the bounded charge will be always 2 regardless of its original value.
	 * @param fragmentationMethodIndex fragmentation method index
	 * @return the bounded charge value
	 */	
	public int getBoundedCharge(int fragmentationMethodIndex) { return charges[fragmentationMethodIndex]; }
	
	/**
	 * Return correct edges with respect to a peptide
	 * @param pep the peptide
	 * @return correct edges
	 */
	public HashSet<Edge> getCorrectEdges(Peptide pep){
		HashSet<Edge> ret = new HashSet<Edge>();
		
		for(ArrayList<Node> k : edgeMap.keySet()){
			Edge e = edgeMap.get(k);
			if(e.isCorrect(pep)) ret.add(e);
		}
		
		return ret;
		
	}
		
	/**
	 * Return correct nodes with respect to a peptide
	 * @param pep the peptide
	 * @return correct nodes
	 */
	public HashSet<Node> getCorrectNodes(Peptide pep){
		HashSet<Node> ret = new HashSet<Node>();
		
		for(Node node : getAllNodes()){
			if(node.isCorrect(pep))
				ret.add(node);
		}
		return ret;
	}
	
	/**
	 * Get an edge between specified nodes. If there is no such an edge, return null.
	 * @param l left node
	 * @param r right node
	 * @return the edge
	 */
	public Edge getEdge(Node l, Node r){
		if(edgeMap == null) return null;
		ArrayList<Node> key = new ArrayList<Node>();
		key.add(l); key.add(r);
		return edgeMap.get(key);
	}
	
	/**
	 * Return filteredNodes 
	 * @return filteredNodes
	 */
	public ArrayList<Node> getFilteredNodes(){return filteredNodes;}
	
	/**
	 * Get linked edges to a node
	 * @param node the node to which edges are linked
	 * @return linked edges
	 */
	public ArrayList<Edge> getLinkedEdges(Node node){
		ArrayList<Edge> edges;
		edges = adjList.get(node);
		
		if(edges == null) return new ArrayList<Edge>();
		else	return edges;
	}
	
	/**
	 * Return nodes within the node tolerance from a mass value
	 * @param mass the mass 
	 * @return nodes within the tolerance
	 */
	public ArrayList<Node> getNodesWithinTolerance(float mass){
		return getNodesWithinTolerance(mass, null);
	}
	
	/**
	 * Return nodes within a tolerance from a mass value
	 * @param mass the mass
	 * @param tol the tolerance 
	 * @return nodes within the tolerance
	 */
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
	
	/**
	 * Return parent mass of this graph
	 * @return parent mass
	 */
	public float getParentMass(){		
		return getSinkNode().getMass() + (float)Composition.H2O;
	}
	
	/**
	 * Return sink node 
	 * @return sink node
	 */
	public Node getSinkNode() {
		return filteredNodes.get(filteredNodes.size()-1);
	}
	
	/**
	 * Return source node 
	 * @return source node
	 */
	public Node getSourceNode() {
		return filteredNodes.get(0);
	}
	
	/**
	 * Merge nodes in spectra of different fragmentation method indices 
	 * @param nodes nodes 
	 */
	private void mergeNodes(ArrayList<Node> nodes){
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
					surv.setLR(surv.getLR() + remov.getLR());					
					nodes.remove(nodes.indexOf(remov));
					i--;
				}	
			}
			prev = surv;
		}
		Collections.sort(nodes);
	}
	
	public String toFileString(){
		String s = "";
		for(Node node : allNodes){
			s += node.toFileString();
		}
		return s;
	}

	/**
	 * Update adjacency list for each node.
	 *
	 * @param enzyme the enzyme
	 * @param updateNodePosteriorAccuracy update node posterior accuracy based on edge accuracies (testing)
	 */	
	private void updateAdjList(Enzyme enzyme, boolean updateNodePosteriorAccuracy){
		adjList = new HashMap<Node, ArrayList<Edge>>();
		edgeMap = new HashMap<ArrayList<Node>, Edge>();	
		adjList.put(filteredNodes.get(0), new ArrayList<Edge>());
		
		for(Node r : filteredNodes){
			ArrayList<Edge> connectingEdges = new ArrayList<Edge>();
			ArrayList<Edge> aaEdges = new ArrayList<Edge>();
			
			for(Node l : filteredNodes){
				if(l.getMass() >= r.getMass()) break;
				if(!adjList.containsKey(l)) continue;
				
				Edge connectingEdge = new Edge(l, r, getParentMass(), aaSet);
			
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
	
				if(connectingEdge.isValid()){
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
			if(!connectingEdges.isEmpty()){				
				adjList.put(r, connectingEdges);
			}
		}		
		if(updateNodePosteriorAccuracy){
			for(int i=1; i<filteredNodes.size()-1;i++){
				filteredNodes.get(i).setPosteriorLRnAccuracy();
			}
		}
		for(ArrayList<Node> key : edgeMap.keySet()){
			edgeMap.get(key).setLRandAccuracy();
		}
	}
}
