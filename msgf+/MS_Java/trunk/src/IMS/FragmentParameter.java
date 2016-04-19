package IMS;

import java.util.ArrayList;

import msutil.Peptide;

public class FragmentParameter {
	private int massIndex, locationIndex, flankingResidueIndex;
	
	
	public FragmentParameter(Peptide peptide, int residueNumber){
		massIndex = getMassIndex(peptide);
		locationIndex = getLocationIndex(peptide, residueNumber);
		flankingResidueIndex = getFlankingResidueIndex(peptide, residueNumber);
	}
	
	static public ArrayList<FragmentParameter> getAllFragmentParameters(){
		ArrayList<FragmentParameter> ret = new ArrayList<FragmentParameter>();
		for(int i=1;i<=3;i++){
			for(int j=1;j<=4;j++){
				for(int k=0;k<=15;k++){
					ret.add(new FragmentParameter(i,j,k));
				}
			}
		}
		return ret;
	}
	
	private FragmentParameter(int massIndex, int locationIndex, int flankingResidueIndex){
		this.massIndex = massIndex;
		this.locationIndex = locationIndex;
		this.flankingResidueIndex = flankingResidueIndex;
	}
	
	private int getMassIndex(Peptide peptide){
		int len = peptide.size();
		if(len < 10) return 1;
		if(len < 20) return 2;
		return 3;
	}
	
	private int getLocationIndex(Peptide peptide, int residueNumber){
		return (int)(peptide.getMass(0, residueNumber)/peptide.getMass()*4 + 1);	
	}
	
	private int getResidueIndex(char residue){
		String residues = "ARNDCQEGHILKMFPSTWYV";
		return residues.indexOf(residue);
	}
	
	private int getFlankingResidueIndex(Peptide peptide, int residueNumber){
		char NtermResidue = peptide.get(residueNumber-1).getResidue();
		char CtermResidue = peptide.get(residueNumber).getResidue();
		return residueTable[getResidueIndex(NtermResidue)][getResidueIndex(CtermResidue)];
	}
	
	private int residueTable[][] = {{0,5, 15,  9,  0,  0,  9,  4,  7, 11, 11,  5,  0,  0,  2, 13, 13,  0,  0, 11},
			{  0,  5, 15,  9,  0 , 0,  9,  4,  7, 11, 11,  5,  0,  0,  2, 13, 13,  0,  0, 11},
		{ 14,  5, 15,  9, 14, 14,  9,4,7, 11, 11,5, 14, 14,2, 13, 13, 14, 14, 11},
				{8,5,8,9,8,8,8,4,7,8,8,5,8,8,2,8,8,8,8,8},
			{0,5, 15,9,0,0,9,4,7, 11, 11,5,0,0,2, 13, 13,0,0, 11},
			{0,5, 15,9,0,0,9,4,7, 11, 11,5,0,0,2, 13, 13,0,0, 11},
			{8,5,8,9,8,8,9,4,7,8,8,5,8,8,2,8,8,8,8,8},
			{3,3,3,3,3,3,3,0,3,3,3,3,3,3,0,3,3,3,3,3},
			{6,5,6,6,6,6,6,4,7,6,6,5,6,6,2,6,6,6,6,6},
			 {10,5, 10,9, 10, 10,9,4,7, 11, 10,5, 10, 10,2, 10, 10, 10, 10, 10},
			 {10,5, 10,9, 10, 10,9,4,7, 11, 11,5, 10, 10,2, 10, 10, 10, 10, 10},
				 { 0,5, 15,9,0,0,9,4,7, 11, 11,5,0,0,2, 13, 13,0,0, 11},
			{0,5, 15,9,0,0,9,4,7, 11, 11,5,0,0,2, 13, 13,0,0, 11},
			{0,5, 15,9,0,0,9,4,7, 11, 11,5,0,0,2, 13, 13,0,0, 11},
			{1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1},
			 {12,5, 12,9, 12, 12,9,4,7, 11, 11,5, 12, 12,2, 13, 12, 12, 12, 11},
			 {12,5, 12,9, 12, 12,9,4,7, 11, 11,5, 12, 12,2, 13, 13, 12, 12, 11},
				 { 0,5, 15,9,0,0,9,4,7, 11, 11,5,0,0,2, 13, 13,0,0, 11},
			{0,5, 15,9,0,0,9,4,7, 11, 11,5,0,0,2, 13, 13,0,0, 11},
			 {10,5, 10,9, 10, 10,9,4,7, 11, 11,5, 10, 10,2, 10, 10, 10, 10, 11}};
	
	public boolean equals(Object o){
		if(this == o) return true;
		if(o instanceof FragmentParameter){
			FragmentParameter other = (FragmentParameter)o;
			return this.massIndex == other.massIndex && 
				this.locationIndex == other.locationIndex &&
				this.flankingResidueIndex == other.flankingResidueIndex;
		}
		return false;
	}
	
	public int hashCode(){
		return new Integer(massIndex).hashCode() * new Integer(locationIndex << 4).hashCode() * new Integer(flankingResidueIndex).hashCode();
	}
	
	public String toFileString(){
		return massIndex+" "+locationIndex+" "+flankingResidueIndex;
	}
	
}
