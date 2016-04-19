package frogDeNovo;

import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Peptide;

public class PrintSeqMass {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//String pep = "HBTPPHGDPGT";
		String pep = "LMDSLKCKISGDC";
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile("/home/kwj/UniNovo/Mods_Cys_H.txt");
		
		//System.out.println(AminoAcidSet.getStandardAminoAcidSet().getAminoAcid('C').getComposition());
		//System.out.println(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys().getAminoAcid('C').getComposition());
		
		System.out.println(new Peptide(pep, aaSet).getMass());
		System.out.println(new Peptide(pep, aaSet).getParentMass()+Composition.H); // M+H

	}

}
