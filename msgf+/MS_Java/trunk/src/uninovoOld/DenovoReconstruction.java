package uninovoOld;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;

import msgf.NominalMass;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.ModifiedAminoAcid;
import msutil.Peptide;
import suffixarray.SuffixArray;


public class DenovoReconstruction {
	private ArrayList<Node> nodes; // include source
	
	private SpectrumGraph graph = null;
	private float accuracy = -1;
	private float LR;
	private String s;
	private int gapNum=-1;
	private int ptmNum = -1;
	//static private float norm = 1;
	//static private SpectrumGraph currentGraph = null;
	
	public DenovoReconstruction(ArrayList<Node> nodes, SpectrumGraph graph){  // nodes include source
		this.graph = graph;
		this.nodes = new ArrayList<Node>();
		for(int i=nodes.size()-1; i>=0; i--){
			this.nodes.add(nodes.get(i));
		}
		calculateAccuracy();
	}
	
	public DenovoReconstruction(ArrayList<Node> nodes){  // nodes include source
		this.nodes = nodes;
	}
	
	public DenovoReconstruction(float Ngap, float Cgap, Peptide annotation, Tolerance tol, Tolerance pmtol){
		this.nodes = new ArrayList<Node>();
		
		nodes.add(new Node(0,tol, annotation.getMass(), 0));
		if(Ngap > 0) nodes.add(new Node(Ngap,tol, annotation.getMass(), 0));
		
		for(float m : annotation.getPRMMasses(true, Ngap)){
			nodes.add(new Node(m,tol,annotation.getMass(),0));
		}
		
		nodes.add(new Node(Ngap + annotation.getMass(),tol,annotation.getMass(),0));
		
		Node sink = new Node(Ngap+Cgap+annotation.getMass(), pmtol, annotation.getMass(), 0);
		sink.setSink(pmtol);
		nodes.add(sink);
	}
	
	DenovoReconstruction(SpectrumGraph graph, ArrayList<Integer> indices){ // masses are sorted in ascending order
		this.graph = graph;
		this.nodes = new ArrayList<Node>();
		
		for(int i : indices){
			this.nodes.add(this.graph.getAllNodes().get(i));
		}
	}
	
	private DenovoReconstruction(DenovoReconstruction one, DenovoReconstruction two, int maxPTMNum){
		this.nodes = new ArrayList<Node>();
		this.graph = one.graph;
		int i=0, j=0;
		
		while(i<one.getNodes().size() || j<two.getNodes().size()){
			Node n1 = one.getNodes().get(i);
			Node n2 = two.getNodes().get(j);
			if(n1.getMass() < n2.getMass()){
				nodes.add(n1);
				i++;
			}else if(n1.getMass() > n2.getMass()){
				nodes.add(n2);
				j++;
			}else{
				nodes.add(n1);
				i++;j++;
			}
		}		
	//	if(this.getPTMNum() > maxPTMNum) this.accuracy = 0f;
	//	else 
			calculateAccuracy();
		//calculateNorm();
		//if(norm > 1)  accuracy /= norm;
	}
	
	public boolean isRedundantWRT(DenovoReconstruction other){
		return this.nodes.containsAll(other.nodes);
	}
	
	public ArrayList<Node> getNodes() {
		return nodes;
	}

	public DenovoReconstruction getUnionDenovoReconstruction(DenovoReconstruction other, int maxPTMNum){
		return new DenovoReconstruction(this, other, maxPTMNum);
	}
	
	public float getAccuracy() { return accuracy;}
	public float getLR() { return LR; }
	
	public boolean isCorrect(Peptide pep){
		if(Math.abs(this.getPeptideMass() - pep.getMass())>2) return false;
		
		boolean isCorrect = true;
		Node prev = nodes.get(0);
		
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Edge edge = null;
			if(graph != null && graph.isEdgeGenerated()){
				edge = graph.getEdge(prev, node);		
				
			}else{
				edge = new Edge(prev, node, nodes.get(nodes.size()-1).getMass(), AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
		//		isCorrect &= prev.isCorrect(pep) && node.isCorrect(pep);
			}
			isCorrect &= edge.isCorrect(pep);
			prev = node;
			if(!isCorrect) break;
		}
		
	
		return isCorrect;
	}
	
	public boolean isCorrectWithMassShiftInGap(Peptide pep, float[] offset){
		
		
		if(this.getPeptideMass() - pep.getMass() < -50 || Math.abs(this.getPeptideMass() - pep.getMass())> 200) return false;
		
		int gapNum = getGapNum();	

		
		for(int g = 0; g<=gapNum; g++){
			boolean isCorrect = true;
			Node prev = nodes.get(0);
			int cn = 0;
			float of = 0;
			//System.out.println("*" + pep + " " + cn + " " + g + "\t" + gapNum);
			for(int i=1; i<nodes.size(); i++){
				Node node = nodes.get(i);
				Edge edge = graph.getEdge(prev, node);
				
				if(edge.isGap()){
					
					cn ++;
					
					isCorrect &= prev.isCorrect(pep, of);
					
					if(cn == g){
						
						of = offset[0] = this.getPeptideMass() - pep.getMass();
						offset[1] = edge.getLeftNode().getMass();
						offset[2] = edge.getRightNode().getMass();
					}
					
					
					isCorrect &= node.isCorrect(pep, of);
					
				}else{
				
					isCorrect &= edge.isCorrect(pep, of);
				}
				
				prev = node;
				
				if(!isCorrect){
					
					break;
				}
			}
			//System.out.println("*" + pep + " " + cn + " " + g + "\t" + gapNum);
			if(isCorrect){
			//	System.out.println("*" + pep + " " + " " + of + "\t" + gapNum);
				return isCorrect;
			}
		}
	
		return false;
		
	}
	
	
	@Override
	public String toString(){
		if(s != null) return s;
		StringBuffer s = new StringBuffer();
		Node prev = nodes.get(0);
		
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Edge edge = null;
			if(graph != null){
				edge = graph.getEdge(prev, node);			
			}else{
				edge = new Edge(prev, node, nodes.get(nodes.size()-1).getMass(), AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
			}
				s.append(edge);
			prev = node;			
		}
		/*
		for(int i=1; i<nodes.size(); i++){
			NewNode node = nodes.get(i);
			s += NominalMass.toNominalMass(node.getMass())+" ";
		}*/
		
		this.s = s.toString();
		return this.s;
	}
	
	public String toFileString(){
		String s = "";
		for(Node node : nodes){
			s += graph.getAllNodes().indexOf(node)+",";
		}
		return s;
	}
	
	
	
	
	private void getUnModifiedGapMassRepresentations(ArrayList<Integer> gaps, int numMod, ArrayList<ArrayList<Integer>> out, AminoAcidSet aaSet){
		if(numMod > aaSet.getMaxNumberOfVariableModificationsPerPeptide()) return;
		if(gaps.size() >= this.nodes.size()-1){
			out.add(gaps);
			return;
		}
		
		Node prev = nodes.get(gaps.size());
		Node node = nodes.get(gaps.size()+1);
		
		Edge edge = null; // TODO opt later...
		if(graph != null){
			edge = graph.getEdge(prev, node);			
		}else{
			edge = new Edge(prev, node, nodes.get(nodes.size()-1).getMass(), AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
		}
		
		
		
		ArrayList<ArrayList<AminoAcid>> modList = edge.getModifiedAAList();
		
		if(modList == null || modList.isEmpty() || edge.isPTM()){
			gaps.add(NominalMass.toNominalMass(edge.getMass()));
			getUnModifiedGapMassRepresentations(gaps, numMod, out, aaSet);
		}else{
			for(ArrayList<AminoAcid> mods : modList){
				ArrayList<Integer> newGaps = new ArrayList<Integer> (gaps);
				int m = NominalMass.toNominalMass(edge.getMass());
				
				float mm = 0;
				for(AminoAcid aa : mods){
					mm += ((ModifiedAminoAcid)aa).getModification().getMass();
				}
				newGaps.add(m - NominalMass.toNominalMass(mm));
				getUnModifiedGapMassRepresentations(newGaps, numMod+mods.size(), out, aaSet);
			
			}
		}
		

	}
	
	
	public ArrayList<ArrayList<Integer>> getGapMassRepresentation(AminoAcidSet aaSet){
		ArrayList<ArrayList<Integer>> out = new ArrayList<ArrayList<Integer>>();
		 getUnModifiedGapMassRepresentations(new ArrayList<Integer>(), 0, out, aaSet);
		 return out;

		/*ArrayList<Integer> ret = new ArrayList<Integer>();
		Node prev = nodes.get(0);
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			
			Edge edge = null; // TODO opt later...
			if(graph != null){
				edge = graph.getEdge(prev, node);			
			}else{
				edge = new Edge(prev, node, nodes.get(nodes.size()-1).getMass(), aaset);
			}
			
			ret.add(NominalMass.toNominalMass(edge.getMass()));
			prev = node;
		}
		return ret;*/
	}
	
	@Override
	public boolean equals(Object other){
		//return this.nodes.equals(((DenovoReconstruction)other).nodes);
		return this.toString().equals(((DenovoReconstruction)other).toString());
	}
	
	public int length() {
		return nodes.size()-1;
	}
	
	public float getPeptideMass() { return (this.nodes.get(nodes.size()-1).getMass());}
	
	public int getPTMNum(){
		if(ptmNum>=0) return ptmNum;
		ptmNum = 0;
		Node prev = nodes.get(0);
		
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Edge edge = graph.getEdge(prev, node);	
			if(edge.isPTM()) ptmNum++;
			prev = node;
		}
		return ptmNum;
	}
	
	public int getGapNum(){
		if(gapNum>=0) return gapNum;
		gapNum = 0;
		Node prev = nodes.get(0);
		
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Edge edge = graph.getEdge(prev, node);	
			if(edge.isGap()) gapNum++;
			prev = node;
		}
		return gapNum;
	}
	
	public float getMass(){
		return graph.getSinkNode().getMass();
	}
	
	private void calculateAccuracy(){
		
		float acc = 0f; // source
		float minEdgeAcc = 1;
		
		LR = nodes.get(0).getLR();
	
		Node prev = nodes.get(0);
		//int len = 0;
		//System.out.println(LR+"*");
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Edge edge = graph.getEdge(prev, node);		
			if(edge == null || !edge.isValid()){
				this.accuracy = 0.001f;
				this.LR = 0f;
				return;
			}
			
			minEdgeAcc = Math.min(minEdgeAcc, edge.getAccuracy());
			LR += edge.getLR();//Edge.getEdge(prev, node).getLR();//node.getLR();//Edge.getEdge(prev, node).getLR();//node.getLR();//
			float w = node.getAccuracy();
		//	System.out.println(LR+"*");
			if(i<nodes.size()-1){
				acc += w;
			//	acc -= + 0.05f * Math.sqrt(w * (1-w));
			}

			prev = node;
		}
		
	//	len = 0;

		for(Edge e : getmST()){
			float w = e.weightForMST();
			acc -= w;
			minEdgeAcc = Math.min(minEdgeAcc, e.getAccuracy());
		//	acc += + 0.005f * Math.sqrt(w * (1-w));
		}
		
		//acc *= 8f/7f;
		//acc *= graph.getSinkNode().getAccuracy();
		
		//if(graph.getSinkNode().getTol().getToleranceAsDa(getMass()) > 0.5f) acc *=0.8f;
		
		acc = Math.max(0.001f, acc);//one.getAccuracy()+two.getAccuracy()-1
		acc = Math.min(1f, acc);
		acc = Math.min(minEdgeAcc, acc);
		
		if(this.length() == 1) this.accuracy = 1;
		else this.accuracy = acc;
	}
	
	/*
	private void calculateNorm(){
		if(currentGraph == graph) return;
		
		currentGraph = graph;
		float acc = 0f; // source

		
		for(int i=1; i<graph.getNodes().size()-1; i++){
			Node node = graph.getNodes().get(i);
			float w = node.getAccuracy();

			acc += w;			


		}

		for(Edge e : getMST()){
			float w = e.getAccuracy();
			acc -= w;
		}
	
		norm = acc;
		
	}
*/
	private ArrayList<Edge> getmST(){
		PriorityQueue<Node> queue = new PriorityQueue<Node>(nodes.size(), new AttachmentComparator());
		ArrayList<Edge> mST = new ArrayList<Edge>();
		ArrayList<Node> nodesNotInTree = new ArrayList<Node>();
		HashMap<Node, Node> connectedNode = new HashMap<Node, Node>();
	
		for(int i=1; i<nodes.size()-1;i++){// init
			Node node = nodes.get(i);
			if(i==1){
				node.setAttachmentCost(0f);
			}else{
				node.setAttachmentCost(Float.MAX_VALUE);
			}
			nodesNotInTree.add(node);
			queue.add(node);
		}
		
		while(!nodesNotInTree.isEmpty()){
			Node node = queue.poll();
			Node cnode = connectedNode.get(node);
			
			if(cnode != null){
				Edge e;
				if(node.getMass() < cnode.getMass()){
					e = graph.getEdge(node, cnode);
				}else{
					e = graph.getEdge(cnode, node);
				}
				mST.add(e);
			}
			nodesNotInTree.remove(node);
			
			for(Node n : nodesNotInTree){
				Edge e;
				if(n.getMass() < node.getMass()){
					e = graph.getEdge(n, node);
				}else{
					e = graph.getEdge(node, n);
				}
				if(e ==null || !e.isValid()) continue;
				float weight = e.weightForMST();
				
				if(weight < n.getAttachmentCost()){// change key
					queue.remove(n);
					n.setAttachmentCost(weight);
					queue.add(n);
					connectedNode.put(n, node);
				}
			}
		}
		return mST;
	}
	
	/*
	private ArrayList<Edge> getMST(){
		PriorityQueue<Node> queue = new PriorityQueue<Node>(graph.getNodes().size(), Collections.reverseOrder( new AttachmentComparator()));
		
		ArrayList<Edge> MST = new ArrayList<Edge>();
		ArrayList<Node> nodesNotInTree = new ArrayList<Node>();
		HashMap<Node, Node> connectedNode = new HashMap<Node, Node>();
	
		for(int i=1; i<graph.getNodes().size()-1;i++){// init
			Node node = graph.getNodes().get(i);
			if(i==1){
				node.setAttachmentCost(0f);
			}else{
				node.setAttachmentCost(-Float.MAX_VALUE);
			}
			nodesNotInTree.add(node);
			queue.add(node);
		}
		
		while(!nodesNotInTree.isEmpty()){
			Node node = queue.poll();
			Node cnode = connectedNode.get(node);
			
			if(cnode != null){
				Edge e;
				if(node.getMass() < cnode.getMass()){
					e = graph.getEdge(node, cnode);
				}else{
					e = graph.getEdge(cnode, node);
				}
				MST.add(e);
			}
			nodesNotInTree.remove(node);
			
			for(Node n : nodesNotInTree){
				Edge e;
				if(n.getMass() < node.getMass()){
					e = graph.getEdge(n, node);
				}else{
					e = graph.getEdge(node, n);
				}
				if(!e.isValid()) continue;
				float weight = e.getAccuracy();
				
				if(weight > n.getAttachmentCost()){// change key
					queue.remove(n);
					n.setAttachmentCost(weight);
					queue.add(n);
					connectedNode.put(n, node);
				}
			}
		}
		return MST;
	}
	*/
	
	public ArrayList<String> getMatches(SuffixArray sa, AminoAcidSet aaSet, boolean allowPTM){
		ArrayList<String> matches = new ArrayList<String>();
		
		String longestAASeq = "";
		String tmp = "";
		
		float p =0, s = 0;
		float pt =0, st = 0;
		
		Node prev = nodes.get(0);
		
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Edge edge = graph.getEdge(prev, node);	
			
			if(!edge.isGap()){
				tmp += edge;
				st = this.getPeptideMass() -  edge.getRightNode().getMass();
			}else{
				tmp = "";
				pt = edge.getRightNode().getMass();
			}
			
			if(longestAASeq.length() <= tmp.length()){
				longestAASeq = tmp.replace('Q', 'K').replace('L', 'I');
				p = pt; s = st;
			}
			
			prev = node;
		}
		
		if(longestAASeq.length() < 3) return matches;
		
		
		ArrayList<String> ms = sa.getAllMatchedStrings(longestAASeq, (int)(p/57), (int)(s/57));//, 
		
		//System.out.println(ms.size());
		
		for(String pep : ms){
			for(int i=0; i<pep.length(); i++){
				String pe = "";
				for(int j=i; j<pep.length(); j++){
					char a = pep.charAt(j);
					if(aaSet.contains(a)) pe += a;
					
						
					if(pe.length()<5 || !pe.replace('Q', 'K').replace('L', 'I').contains(longestAASeq)) continue; 
					
					String toadd = "";
					if(allowPTM){						
						float[] offset = new float[3];
						if(isCorrectWithMassShiftInGap(new Peptide(pe, aaSet), offset)){									
							toadd = pe + "\t" + offset[0] + "\t" + offset[1] + "\t" + offset[2];
						}
					}else{
						if(isCorrect(new Peptide(pe, aaSet))){
							toadd = pe;
						}
					}
					if(!toadd.isEmpty() && !matches.contains(toadd)) matches.add(toadd);
				}
			}
			
			
			//choose min offset one?
			
			
			
		}
		
		/*
		float minoff = Float.MAX_VALUE;
		for(int i=0; i<matches.size(); i++){
			String[] token = matches.get(i).split("\t");
			if(token.length < 2) continue;
			minoff = Math.min(minoff, Math.abs(Float.parseFloat(token[1])));
		}
		
		for(int i=0; i<matches.size(); i++){
			String[] token = matches.get(i).split("\t");
			if(token.length < 2) continue;
			if(Math.abs(Float.parseFloat(token[1])) > minoff)
				matches.remove(i--);
		}
		*/
		
		return matches;
		
	}
	
	
	private class AttachmentComparator implements Comparator<Node>{
		@Override
		public int compare(Node o1, Node o2) {
			return new Float(o1.getAttachmentCost()).compareTo(new Float(o2.getAttachmentCost()));
		}
	}
	
	public static class DenovoReconstructionAccuracyComparator implements Comparator<DenovoReconstruction>{ 
		@Override
		public int compare(DenovoReconstruction arg0, DenovoReconstruction arg1) {
			int r = new Float(arg0.accuracy).compareTo(new Float(arg1.accuracy));
			if(r != 0) return r;
			
			return -new Integer(arg0.length()).compareTo(new Integer(arg1.length()));
		}
		
		static public DenovoReconstructionAccuracyComparator get(){
			return new DenovoReconstructionAccuracyComparator();
		}
	}
	
	public static class DenovoReconstructionLRComparator implements Comparator<DenovoReconstruction>{ // 
		@Override
		public int compare(DenovoReconstruction arg0, DenovoReconstruction arg1) {
			return new Float(arg0.LR).compareTo(new Float(arg1.LR));
		}
		
		static public DenovoReconstructionLRComparator get(){
			return new DenovoReconstructionLRComparator();
		}
	}
}
