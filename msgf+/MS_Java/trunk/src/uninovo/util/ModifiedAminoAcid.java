package uninovo.util;

// for variable modification
public class ModifiedAminoAcid extends AminoAcid {
	private boolean isResidueSpecific;
	private Modification mod;
	private char unmodResidue;
	
	public ModifiedAminoAcid(AminoAcid aa, Modification mod, char residue)
	{
		super(residue, mod.getName()+" "+aa.getName(), aa.getAccurateMass()+mod.getAccurateMass());
		this.mod = mod;
		this.unmodResidue = aa.getResidue();
		super.setProbability(aa.getProbability());
		if(Character.isUpperCase(unmodResidue))
			isResidueSpecific = true;
	}
	
	public Modification getModification()	{ return mod; }
	@Override
	public String getResidueStr()
	{
		StringBuffer buf = new StringBuffer();
		buf.append(unmodResidue);
		float modMass = mod.getMass();
		if(modMass >= 0)
			buf.append('+');
		buf.append(String.format("%.3f", modMass));
		return buf.toString();
	}
	public char getUnmodResidue() 	{ return unmodResidue; }
	
	@Override
	public boolean isModified()		{ return true; }
	
	public boolean isResidueSpecific()	{ return isResidueSpecific; }
}
