package uninovoOld.analysis;

import java.io.IOException;
import java.util.ArrayList;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.WindowFilter;
import unexplainedpeak.GenerateIonDistGraph;
import uninovoOld.IPDGenerator;

public class RankTest {
	static public void main(String[] args) throws IOException{
		WindowFilter filter = new WindowFilter(6, 50);
		int charge =-1; 
		int emphasizeIthIon = -1;
		//int maxCharge = 2;
	//	Tolerance tol = new Tolerance(20f, true);
		Tolerance pmtol = new Tolerance(20f, true);
		Tolerance tol = new Tolerance(0.5f, false);
		//Tolerance pmtol = new Tolerance(30f, true);
		//Tolerance tol = new Tolerance(100f, true);
	//"/home/kwj/workspace/inputs/Training/HCD_train.mgf";
		String para = "/home/kwj/Dropbox/PAR/ETDTrypsin_train.parsnew";//"/home/kwj/workspace/inputs/nrp.mgf";//
		String inputmgf = System.getProperty("user.home") + "/Dropbox/Test/ETDTrypsin_test.mgf";// trainmgf;//Zubarev_HCD_Annotated.mgf
		//String inputmgf =System.getProperty("user.home") + "/workspace/inputs/Zubarev_HCD_Annotated.mgf";// trainmgf;//Zubarev_HCD_Annotated.mgf
		
		String outmgf = inputmgf.substring(0, inputmgf.lastIndexOf(".")) + "_ranked_" + emphasizeIthIon + ".mgf";
	//	String pepoutmgf = "/home/kwj/workspace/inputs/Training/ETD_Tryp_Confident_PepNovo+"+charge+".mgf";

		int maxIterationNum = 100;

		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		IPDGenerator gen = new IPDGenerator(para, aaSet, 0).filter(filter);
		gen.run(inputmgf, outmgf, charge, maxIterationNum, 0, false, emphasizeIthIon);
		
		//charge = 3;

		ArrayList<IonType> sigIons= gen.getSigIonsOrderedByIntensityWithOutNoiseIon(charge);
		sigIons.remove(sigIons.size()-1);
		float probabiltiy = 0.1f;
	
		sigIons = GenerateIonDistGraph.gen(inputmgf, 150, sigIons, charge, tol, probabiltiy);

		sigIons = GenerateIonDistGraph.gen(outmgf, 150, sigIons, charge, tol, probabiltiy);
		

	}
}
