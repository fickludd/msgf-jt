package uninovo.util;

import java.util.ArrayList;


public interface DeNovoNodeFactory<T extends Matter> {
	public boolean contains(T node);
	public AminoAcidSet getAASet();
	public T getComplementNode(T srm, T pmNode);
	public DeNovoGraph.Edge<T> getEdge(T curNode, T prevNode);
	public ArrayList<DeNovoGraph.Edge<T>> getEdges(T curNode);
	public Enzyme getEnzyme();
	public T getInfinity();
	public ArrayList<T> getIntermediateNodes(ArrayList<T> destNodes);
	public T getNextNode(T curNode, AminoAcid aa);
	public ArrayList<T> getNodes(float peptideMass, Tolerance tolerance);
	public T getPreviousNode(T curNode, AminoAcid aa);
	public T getZero();
	public boolean isReverse();
	public int size();
	public Sequence<T> toCumulativeSequence(boolean isPrefix, Peptide pep);
}
