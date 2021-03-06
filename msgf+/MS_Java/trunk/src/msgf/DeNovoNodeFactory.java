package msgf;

import java.util.ArrayList;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Matter;
import msutil.Peptide;
import msutil.Sequence;

public interface DeNovoNodeFactory<T extends Matter> {
	public AminoAcidSet getAASet();
	public T getZero();
	public T getInfinity();
	public ArrayList<T> getNodes(float peptideMass, Tolerance tolerance);
	public T getComplementNode(T srm, T pmNode);
	public ArrayList<T> getIntermediateNodes(ArrayList<T> destNodes);
	public ArrayList<DeNovoGraph.Edge<T>> getEdges(T curNode);
	public DeNovoGraph.Edge<T> getEdge(T curNode, T prevNode);
	public Sequence<T> toCumulativeSequence(boolean isPrefix, Peptide pep);
	public T getPreviousNode(T curNode, AminoAcid aa);
	public T getNextNode(T curNode, AminoAcid aa);
	public int size();
	public boolean contains(T node);
	public boolean isReverse();
	public Enzyme getEnzyme();
}
