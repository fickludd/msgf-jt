package uninovo.util;

import java.util.ArrayList;


public abstract class DeNovoGraph<T extends Matter> {
	public static class Edge<T extends Matter> {
		// scores
		private int cleavageScore;
		private int errorScore;
		private int index;	
		private float mass;

		private T prevNode;
		private float probability;
		
		public Edge(T prevNode, float probability, int index, float mass) 
		{
			this.prevNode = prevNode;
			this.probability = probability;
			this.index = index;
			this.mass = mass;
		}
		public int getEdgeIndex() {
			return index;
		}
		public float getEdgeMass() {
			return mass;
		}
		public float getEdgeProbability() {
			return probability;
		}
		public int getEdgeScore() {
			return cleavageScore+errorScore;
		}
		public T getPrevNode() {
			return prevNode;
		}
		public void setCleavageScore(int cleavageScore)
		{
			this.cleavageScore = cleavageScore;
		}
		public void setEdgeMass(float mass)
		{
			this.mass = mass;
		}
		public void setErrorScore(int errorScore)
		{
			this.errorScore = errorScore;
		}
	}
	protected T destination;
	protected ArrayList<T> intermediateNodes;
	protected T pmNode;
	protected ArrayList<T> sinkNodes;
	
	protected T source;
	public abstract AminoAcidSet getAASet();
	public abstract T getComplementNode(T node);
	public T getDestination() { return destination; }
	//	public abstract int getEdgeScore(T curNode, T prevNode);
	public abstract ArrayList<Edge<T>> getEdges(T curNode);
	
	public ArrayList<T> getIntermediateNodeList() { return intermediateNodes; }
	public abstract int getNodeScore(T node);
	public T getPMNode() { return pmNode; }
	public abstract int getScore(Annotation annotation);
public abstract int getScore(Peptide pep);
	public ArrayList<T> getSinkList() { return sinkNodes; }
	public T getSource() { return source; }
	
	public abstract boolean isReverse();
}
