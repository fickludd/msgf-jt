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
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;

import uninovo.util.AminoAcidSet;

/**
 * The Class DenovoSequence represents de novo sequence
 *
 * @author kyowon Jeong
 */
public class DenovoSequence {
	
	/**
	 * The Class AttachmentComparator compares the attachment cost of nodes
	 *
	 * @author Kyowon Jeong
	 */
	private class AttachmentComparator implements Comparator<Node>{
		@Override
		public int compare(Node o1, Node o2) {
			return new Float(o1.getAttachmentCost()).compareTo(new Float(o2.getAttachmentCost()));
		}
	}
	
	/**
	 * The Class DenovoSequenceLRComparator compares the LR scores of de novo sequences
	 *
	 * @author Kyowon Jeong
	 */
	public static class DenovoSequenceLRComparator implements Comparator<DenovoSequence>{ // 
		/**
		 * Gets the.
		 *
		 * @return the de novo sequence lr comparator
		 */
		static public DenovoSequenceLRComparator get(){
			return new DenovoSequenceLRComparator();
		}
		
		@Override
		public int compare(DenovoSequence arg0, DenovoSequence arg1) {
			return new Float(arg0.LR).compareTo(new Float(arg1.LR));
		}
	}
	
	/** The accuracy of this. */
	private float accuracy = -1;
	
	/** The gap num. */
	private int gapNum=-1;
	
	/** The spectrum graph on which sequences are generated. */
	private SpectrumGraph graph = null;
	
	/** The lr score. */
	private float LR;
	
	/** The nodes. */
	private ArrayList<Node> nodes; // include source
	
	/** The number of PTMs for future use. */
	private int ptmNum = -1;
	
	/** The string representation of this. */
	private String stringRepresentation;
	
	/**
	 * Instantiates a new de novo sequence. Accuracy is computed.
	 *
	 * @param nodes the nodes
	 * @param graph the graph
	 */
	public DenovoSequence(ArrayList<Node> nodes, SpectrumGraph graph){  // nodes include source
		this.graph = graph;
		this.nodes = new ArrayList<Node>();
		for(int i=nodes.size()-1; i>=0; i--){
			this.nodes.add(nodes.get(i));
		}
		calculateAccuracy();
	}
	
	/**
	 * Instantiates a new denovo sequence by taking union of nodes in two sequences.
	 *
	 * @param one the first seq
	 * @param two the second seq
	 */
	private DenovoSequence(DenovoSequence one, DenovoSequence two){
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
		calculateAccuracy();		
	}

	/**
	 * Calculate accuracy.
	 */
	private void calculateAccuracy(){
		
		float acc = 0f; // source
		float minEdgeAcc = 1;
		
		LR = nodes.get(0).getLR();
	
		Node prev = nodes.get(0);
		for(int i=1; i<nodes.size(); i++){
			Node node = nodes.get(i);
			Edge edge = graph.getEdge(prev, node);		
			if(edge == null || !edge.isValid()){
				this.accuracy = 0.001f;
				this.LR = 0f;
				return;
			}
			
			minEdgeAcc = Math.min(minEdgeAcc, edge.getAccuracy());
			LR += edge.getLR();
			float w = node.getAccuracy();
	
			if(i<nodes.size()-1){
				acc += w;	
			}

			prev = node;
		}
		for(Edge e : getmST()){
			float w = e.weightForMST();
			acc -= w;
			minEdgeAcc = Math.min(minEdgeAcc, e.getAccuracy());	
		}
		
		acc = Math.max(0.001f, acc);
		acc = Math.min(1f, acc);
		acc = Math.min(minEdgeAcc, acc);
		
		if(this.length() == 1) this.accuracy = 1;
		else this.accuracy = acc;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object other){
		return this.toString().equals(((DenovoSequence)other).toString());
	}
	
	/**
	 * Gets the accuracy.
	 *
	 * @return the accuracy
	 */
	public float getAccuracy() { return accuracy;}
			
	/**
	 * Gets the gap number.
	 *
	 * @return the gap number
	 */
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
	
	/**
	 * Gets the lr score.
	 *
	 * @return the lr score
	 */
	public float getLR() { return LR; }
	
	/**
	 * Gets the mass.
	 *
	 * @return the mass
	 */
	public float getMass(){
		return graph.getSinkNode().getMass();
	}
	
	/**
	 * Gets the minimum spanning tree (to get the Hunters bound).
	 *
	 * @return the mst
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

	/**
	 * Gets the nodes.
	 *
	 * @return the nodes
	 */
	public ArrayList<Node> getNodes() {
		return nodes;
	}
	
	/**
	 * Gets the PTM number.
	 *
	 * @return the PTM num
	 */
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
	
	/**
	 * Gets the union denovo sequence.
	 *
	 * @param other the other
	 * @return the union denovo sequence
	 */
	public DenovoSequence getUnionDenovoSequence(DenovoSequence other){
		return new DenovoSequence(this, other);
	}
	
	/**
	 * Checks if this is redundant wrt other sequence.
	 *
	 * @param other the other
	 * @return true, if is redundant 
	 */
	public boolean isRedundantWRT(DenovoSequence other){
		return this.nodes.containsAll(other.nodes);
	}
	
	/**
	 * Length.
	 *
	 * @return the length of this
	 */
	public int length() {
		return nodes.size()-1;
	}	
	
	/**
	 * To file string.
	 *
	 * @return the string
	 */
	public String toFileString(){
		String s = "";
		for(Node node : nodes){
			s += graph.getAllNodes().indexOf(node)+",";
		}
		return s;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString(){
		if(stringRepresentation != null) return stringRepresentation;
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
	
		this.stringRepresentation = s.toString();
		return this.stringRepresentation;
	}
}
