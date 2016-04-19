package msscorer;

import msgf.NominalMass;
import msgf.NominalMassFactory;

// Fast scorer for DB search, consider edges
public class DBScanScorer extends FastScorer {

	private float[] nodeMass = null;
	private NewRankScorer scorer = null;
	private Partition partition;
	private float probPeak;
	private boolean mainIonDirection;	// prefix: true, suffix: false
	
	public DBScanScorer(NominalMassFactory factory,	NewScoredSpectrum<NominalMass> scoredSpec, int peptideMass) 
	{
		super(factory, scoredSpec, peptideMass);
		this.scorer = scoredSpec.getScorer();
		
		nodeMass = new float[peptideMass];
		
		for(int i=0; i<nodeMass.length; i++)
			nodeMass[i] = -1;
		
		// assign node mass
		nodeMass[0] = 0;
		for(int nominalMass=1; nominalMass<nodeMass.length; nominalMass++)
		{
			NominalMass node = factory.getInstanceOfIndex(nominalMass);
			if(node != null)
				nodeMass[nominalMass] = scoredSpec.getNodeMass(node);
		}
		
		mainIonDirection = scoredSpec.getMainIon().isPrefixIon();
		
		partition = scoredSpec.getPartition();
		probPeak = scoredSpec.getProbPeak();
	}

	// fromIndex: inclusive, toIndex: exclusive
	@Override
	public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr, int fromIndex, int toIndex)
	{
		int nodeScore = super.getScore(prefixMassArr, nominalPrefixMassArr, fromIndex, toIndex);
		int edgeScore = 0;
		if(!mainIonDirection)
		{
			int nominalPeptideMass = nominalPrefixMassArr[toIndex-1];
			for(int i=toIndex-2; i>=fromIndex; i--)
				edgeScore += getEdgeScoreInt(nominalPeptideMass-nominalPrefixMassArr[i], nominalPeptideMass-nominalPrefixMassArr[i+1], (float)(prefixMassArr[i+1]-prefixMassArr[i]));
		}
		else
		{
			for(int i=fromIndex+1; i<=toIndex-1; i++) 
				edgeScore += getEdgeScoreInt(nominalPrefixMassArr[i], nominalPrefixMassArr[i-1], (float)(prefixMassArr[i]-prefixMassArr[i-1]));
		}
		
		return nodeScore + edgeScore;
	}
	
	@Override
	public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
		return getEdgeScoreInt(curNode.getNominalMass(), prevNode.getNominalMass(), theoMass);
	}
	
	private int getEdgeScoreInt(int curNominalMass, int prevNominalMass, float theoMass) {
		int ionExistenceIndex = 0;
		float curMass = nodeMass[curNominalMass];
		if(curMass >= 0)
			ionExistenceIndex += 1;
		float prevMass = nodeMass[prevNominalMass];
		if(prevMass >= 0)
			ionExistenceIndex += 2;
		
		float edgeScore = scorer.getIonExistenceScore(partition, ionExistenceIndex, probPeak);
		if(ionExistenceIndex == 3)
			edgeScore += scorer.getErrorScore(partition, curMass-prevMass-theoMass);
		return Math.round(edgeScore);
	}
}
