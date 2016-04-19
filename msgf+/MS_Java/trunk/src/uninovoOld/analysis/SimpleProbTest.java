package uninovoOld.analysis;

import java.io.IOException;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;
import uninovoOld.train.UniNovoTrainer;

public class SimpleProbTest {

	static public void main(String[] args) throws IOException{
			
			int charge = 3; 
			Tolerance tol = new Tolerance(0.5f, false);
			Tolerance pmtol = new Tolerance(2.0f, false);
			String trainmgf = "/home/kwj/workspace/inputs/Training/CID_Tryp_Confident.mgf";
			String para = trainmgf + charge + ".parSimple";
		
			int maxIonRank = 40;
			int maxRank = 10;
			int maxIterationNum = 1;
			
			AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
			PeakParameter.setGroupNum(1);
			PeakParameter.setPartitionNum(1);
			SpectrumParameter.setSpecMzRangeNum(1);
			PeakParameter.setPeakIntensityRatioNum(1);
			
			UniNovoTrainer trainer = new UniNovoTrainer(trainmgf, tol, pmtol, aaSet).doNotDiscard();
			trainer.train(para, charge, maxIterationNum, Enzyme.TRYPSIN);
	
		}
}
