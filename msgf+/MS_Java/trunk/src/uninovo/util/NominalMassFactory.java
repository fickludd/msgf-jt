package uninovo.util;

import java.util.ArrayList;


public class NominalMassFactory extends MassFactory<NominalMass> {
	//	// static methods
	private static NominalMassFactory defaultNominalMassFactory = new NominalMassFactory(50);
	public static NominalMass getInstanceFor(float mass)
	{
		return defaultNominalMassFactory.getInstance(mass);
	}
	private int[] aaMassIndex;
	private NominalMass[] factory;
	private NominalMass inf;
	
	private float rescalingConstant = Constants.INTEGER_MASS_SCALER;

	private NominalMass zero;
	
	public NominalMassFactory(AminoAcidSet aaSet, Enzyme enzyme, int maxLength)
	{
		super(aaSet, enzyme, maxLength);
		int heaviestNominalMass = aaSet.getHeaviestAA().getNominalMass();
		int maxIndex = heaviestNominalMass*maxLength;
		factory = new NominalMass[maxIndex+2];
		zero = factory[0] = new NominalMass(0);
		inf = factory[factory.length-1] = new NominalMass(factory.length-1);
		aaMassIndex = new int[128];
		for(AminoAcid aa : aaSet)
			aaMassIndex[aa.getResidue()] = aa.getNominalMass();
		makeAllPossibleMasses(true);
	}

	private NominalMassFactory(int maxLength)
	{
		super(null, null, maxLength);
	}

	@Override
	public boolean contains(NominalMass node) {
		int index = node.getNominalMass();
		if(index < 0 || index >= factory.length)
			return false;
		return factory[index] != null;
	}
	
	@Override
	public NominalMass getComplementNode(NominalMass srm, NominalMass pmNode) {
		int index = pmNode.getNominalMass() - srm.getNominalMass();
		if(factory[index] != null)
			return factory[index];
		else
			return new NominalMass(index);
	}
	
	@Override
	public ArrayList<DeNovoGraph.Edge<NominalMass>> getEdges(NominalMass curNode)
	{
		return edgeMap.get(curNode);
//		if(edgeMap != null)
//			return edgeMap.get(curNode);
//		int curIndex = curNode.getNominalMass();
//		ArrayList<DeNovoGraph.Edge<NominalMass>> edges = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
//		for(AminoAcid aa : aaSet)
//		{
//			int prevIndex = curIndex - aaMassIndex[aa.getResidue()];
//			DeNovoGraph.Edge<NominalMass> edge = new DeNovoGraph.Edge<NominalMass>(new NominalMass(prevIndex), aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
//			int cleavageScore = 0;
//			if(prevIndex == 0 && enzyme != null)
//			{
//				if(enzyme.isCleavable(aa))
//					cleavageScore += enzyme.getPeptideCleavageCredit();
//				else
//					cleavageScore += enzyme.getPeptideCleavagePenalty();
//			}
//			edge.setCleavageScore(cleavageScore);
//			edges.add(edge);
//		}
//		return edges;
	}
	
	@Override
	public NominalMass getInfinity() {
		return inf;
	}
	
	public NominalMass getInstance(float mass)
	{
		int massIndex = getMassIndex(mass);
		return getInstanceOfIndex(massIndex);
	}
	
	// returns instance exists in the factory
	public NominalMass getInstanceOfIndex(int index)
	{
		if(index < factory.length)
			return factory[index];
		else
			return null;
	}

	public float getMassFromIndex(int massIndex)
	{
		return massIndex/rescalingConstant;
	}

	public int getMassIndex(float mass)
	{
		return Math.round(mass*rescalingConstant);
	}
	
	@Override
	public NominalMass getNextNode(NominalMass curNode, AminoAcid aa) {
		int index = curNode.getNominalMass() + aa.getNominalMass();
		if(factory[index] == null)
			factory[index] = new NominalMass(index);
		return factory[index];
	}
	
	@Override
	public ArrayList<NominalMass> getNodes(float peptideMass, Tolerance tolerance) {
		ArrayList<NominalMass> nodes = new ArrayList<NominalMass>();
		float tolDa = tolerance.getToleranceAsDa(peptideMass);
		int minIndex = getMassIndex(peptideMass-tolDa);
		int maxIndex = getMassIndex(peptideMass+tolDa);
		for(int index = minIndex; index<=maxIndex; index++)
		{
			if(factory[index] != null)
				nodes.add(factory[index]);
		}
		return nodes;
	}

	@Override
	public NominalMass getPreviousNode(NominalMass curNode, AminoAcid aa) 
	{
		int index = curNode.getNominalMass() - aa.getNominalMass();
		if(index < 0)
			return null;
//		if(factory[index] == null)
//			factory[index] = new NominalMass(index);
		return factory[index];
	}
	
public float getRescalingConstant()	{ return rescalingConstant; }
	@Override
	public NominalMass getZero() {
		return zero;
	}
}

