package uninovo.util;

import java.util.ArrayList;
import java.util.Collections;


/**
 * A factory class instantiate compositions.
 * @author sangtaekim
 *
 */
public class CompositionFactory extends MassFactory<Composition>{
	
	private static final int arraySize = 1 << 27;
	private static final int indexMask = 0xFFFFFFE0;
	private static final int offsetMask = 0x0000001F;
	
	public static void main(String[] argv)
	{
	}
	private int[] data;
	private int[] map;
	
	private ArrayList<Composition> tempData;	// temporary
	
	public CompositionFactory(AminoAcidSet aaSet, Enzyme enzyme, int maxLength)
	{
		super(aaSet, enzyme, maxLength);
		this.map = new int[arraySize];
		tempData = new ArrayList<Composition>();
		makeAllPossibleMasses();
	}
	
	// private class for getIntermediateCompositions, don't generate all possible nodes
	private CompositionFactory(AminoAcidSet aaSet, int maxLength)
	{
		super(aaSet, null, maxLength);
		this.map = new int[arraySize];
		tempData = new ArrayList<Composition>();
	}

	protected void add(int number)
	{
		tempData.add(new Composition(number)); 
	}
	
	protected void clear(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		map[index] &= ~(1 << offset);
	}

	@Override
	public boolean contains(Composition node) {
		return isSet(node.number);
	}
	
	private CompositionFactory finalizeCompositionSet()	
	{ 
		if(tempData != null) Collections.sort(tempData);
		data = new int[tempData.size()];
		for(int i=0; i<tempData.size(); i++)
			data[i] = tempData.get(i).getNumber();
		tempData = null;
		return this;
	}
	
	@Override
	public Composition getComplementNode(Composition srm, Composition pmNode) {
		return pmNode.getSubtraction(srm);
	}
	
	public int[] getData()	{ return data; }
	
	@Override
	public ArrayList<DeNovoGraph.Edge<Composition>> getEdges(Composition curNode) {
		// prevNode, score, prob, index
		int curNum = curNode.number;
		ArrayList<DeNovoGraph.Edge<Composition>> edges = new ArrayList<DeNovoGraph.Edge<Composition>>();
		for(AminoAcid aa : aaSet)
		{
			int prevNum = curNum - aa.getComposition().number;
			DeNovoGraph.Edge<Composition> edge = new DeNovoGraph.Edge<Composition>(new Composition(prevNum), aa.getProbability(), aaSet.getIndex(aa), aa.getMass()); 
			if(prevNum == 0 && enzyme != null)
			{
				if(enzyme.isCleavable(aa))
					edge.setCleavageScore(enzyme.getPeptideCleavageCredit());
				else
					edge.setCleavageScore(enzyme.getPeptideCleavagePenalty());
			}
			edges.add(edge);
		}
		return edges;
	}

	@Override
	public Composition getInfinity() {
		return Composition.INF;
	}
	
	// return set of compositions contained in paths from (0,0,0,0,0) to despCompositions
	public ArrayList<Composition> getIntermediateCompositions(Composition source, ArrayList<Composition> destCompositionList)
	{
		CompositionFactory intermediateCompositions = new CompositionFactory(this.aaSet, maxLength);
		
		for(Composition c : destCompositionList)
		{
			intermediateCompositions.setAndAddIfNotExist(c.number);
		}
		
		int start = 0;
		while(true)
		{
			int end = intermediateCompositions.size();
			for(int i=start; i<end; i++)
			{
				int number = intermediateCompositions.tempData.get(i).getNumber();
				for(AminoAcid aa : aaSet)
				{
					Composition aaComp = aa.getComposition(); 
					int prevNumber = number-aaComp.getNumber();
					if(this.isSet(prevNumber) && !intermediateCompositions.isSet(prevNumber))
					{
						intermediateCompositions.setAndAddIfNotExist(prevNumber);
					}
				}
			}
			if(end == intermediateCompositions.size())
				break;
			start = end;
		}
		
		Collections.sort(intermediateCompositions.tempData);
		return intermediateCompositions.tempData;
	}
	
	@Override
	public ArrayList<Composition> getIntermediateNodes(ArrayList<Composition> destCompositionList)
	{
		return getIntermediateCompositions(new Composition(0), destCompositionList);
	}
	
	@Override
	public Composition getNextNode(Composition curNode, AminoAcid aa) {
		int num = curNode.number + aa.getComposition().number;
		return new Composition(num);
	}
	
	@Override
	public ArrayList<Composition> getNodes(float mass, Tolerance tolerance)
	{
		ArrayList<Composition> compositions = new ArrayList<Composition>();
	
		float toleranceDa = tolerance.getToleranceAsDa(mass);
		float minMass = mass-toleranceDa;
		float maxMass = mass+toleranceDa;
		// binary search
		int minIndex=0, maxIndex=data.length, i=-1;
		while(true)
		{
			i = (minIndex+maxIndex)/2; 
			double m = Composition.getMonoMass(data[i]);
			if(m < minMass)
				minIndex = i;
			else if(m > maxMass)
				maxIndex = i;
			else
				break;
			if(maxIndex - minIndex <= 1)
				break;
		}
		for(int cur=i; cur>=0; cur--)
		{
			double m = Composition.getMonoMass(data[cur]);
			if(m >= minMass && m <= maxMass)
				compositions.add(new Composition(data[cur]));
			else if(m < minMass)
				break;
		}
		for(int cur=i+1; cur<data.length; cur++)
		{
			double m = Composition.getMonoMass(data[cur]);
			if(m >= minMass && m <= maxMass)
				compositions.add(new Composition(data[cur]));
			else if(m > maxMass)
				break;
		}
		Collections.sort(compositions);
		return compositions;
	}
	
	@Override
	public Composition getZero() {
		return Composition.NIL;
	}
	
	private boolean isSet(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		if((map[index] & (1 << offset)) == 0)
			return false;
		else
			return true;
	}
	
	protected void makeAllPossibleMasses()
	{
		setAndAddIfNotExist(0);

		Composition[] aaComposition = new Composition[aaSet.size()];
		int index=0;
		for(AminoAcid aa : aaSet)
			aaComposition[index++] = aa.getComposition();

		int start = 0;
		for(int l=0; l<maxLength; l++)
		{
			int end = tempData.size();
			for(int i=start; i<end; i++)
			{
				for(int j=0; j<aaComposition.length; j++)
					setAndAddIfNotExist(tempData.get(i).getNumber()+aaComposition[j].getNumber());
			}
			start = end;
		}
		finalizeCompositionSet();
	}
	
	protected void set(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		map[index] |= (1 << offset);
	}
	

	private void setAndAddIfNotExist(int number)
	{
		int index = (number & indexMask) >>> 5;
		int offset = number & offsetMask;
		if((map[index] & (1 << offset)) == 0)	// nonexistant 
		{
			map[index] |= (1 << offset);	// set
			tempData.add(new Composition(number));	// add
		}
	}	

	@Override
	public int size()	
	{ 
		if(data == null)	// not finalized yet
		{
			if(tempData == null)
				return -1;
			else
				return tempData.size();
		}
		else 
			return data.length; 
	}


}