package uninovo.util;

public class Annotation {
	private AminoAcid nextAA;
	private Peptide peptide;
	private AminoAcid prevAA;
	
	public Annotation(AminoAcid prevAA, Peptide peptide, AminoAcid nextAA)
	{
		this.prevAA = prevAA;
		this.peptide = peptide;
		this.nextAA = nextAA;
	}
	
	public Annotation(String annotationStr, AminoAcidSet aaSet)
	{
		String pepStr = annotationStr.substring(annotationStr.indexOf('.')+1, annotationStr.lastIndexOf('.'));
		char prevAAResidue = annotationStr.charAt(0);
		char nextAAResidue = annotationStr.charAt(annotationStr.length()-1);
		
		prevAA = aaSet.getAminoAcid(prevAAResidue);
		peptide = aaSet.getPeptide(pepStr);
		nextAA = aaSet.getAminoAcid(nextAAResidue);
	}
	
	public AminoAcid getNextAA() {
		return nextAA;
	}

	public Peptide getPeptide() {
		return peptide;
	}

	public AminoAcid getPrevAA() {
		return prevAA;
	}

	public void setNextAA(AminoAcid nextAA) {
		this.nextAA = nextAA;
	}

	public void setPeptide(Peptide peptide) {
		this.peptide = peptide;
	}

	public void setPrevAA(AminoAcid prevAA) {
		this.prevAA = prevAA;
	}
	
	@Override
	public String toString()
	{
		if(peptide == null)
			return null;
		StringBuffer output = new StringBuffer();
		if(prevAA != null)
			output.append(prevAA.getResidueStr());
		output.append("."+peptide.toString()+".");
		if(nextAA != null)
			output.append(nextAA.getResidueStr());
		return output.toString();
	}
}
