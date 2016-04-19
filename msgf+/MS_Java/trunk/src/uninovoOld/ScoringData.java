package uninovoOld;

import java.util.ArrayList;



public class ScoringData{
	private int specID;
	private SpectrumGraph graph;
	private ArrayList<DenovoReconstruction> recs =null;
	private String specfilename;
	//private Tolerance tol, pmtol;
	
	ScoringData(int specID, SpectrumGraph graph, ArrayList<DenovoReconstruction> recs, String specfilename){
		this.specID = specID;
		this.graph = graph;
		this.recs = recs;		
		this.specfilename = specfilename;
	}

	public int getSpecID() {
		return specID;
	}

	public SpectrumGraph getGraph() {
		return graph;
	}

	public String getSpecfilename() {
		return specfilename;
	}

	public ArrayList<DenovoReconstruction> getRecs(){
		return recs;
	}
	/*
	float getScore(Peptide peptide, float offset, int position){ // peptide is not modified. position 0 ~ peptide length -1
		boolean isCorrect = false;
		
		
		
		if(offset !=0){
			String s = peptide.toString().substring(0, position+1) + "+" + offset + peptide.toString().substring(position+1, peptide.size());
			peptide = new Peptide(s, graphs.getAminoAcidSet());
		}else{
			peptide = new Peptide(peptide.toString(), graphs.getAminoAcidSet());
		}

		if(isCorrect) return graphs.getScore(peptide);
		else return 0;
	}*/
	
	
	
	
}
