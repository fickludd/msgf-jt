package uninovoOld.analysis;

import java.util.ArrayList;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.WindowFilter;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;
import uninovoOld.train.FeatureSelector;
import uninovoOld.train.FeatureStatisticsTrainer;
import uninovoOld.train.IonSelector;

public class PrintSigFeatures {
	public static void main(String[] args){
		int charge = 3;
		
		String mgf = "/home/kwj/workspace/inputs/Training/HCD_train.mgf";
		
		if(args.length>1){
			mgf = args[0];
			charge = Integer.parseInt(args[1]);
		}
		
		String para = "/home/kwj/Dropbox/UniNovo/Feature/"+mgf.substring(mgf.lastIndexOf('/'))+ charge + "_.paralatex";
		
		Tolerance tol = new Tolerance(0.5f, false), 
		pmtol = new Tolerance(20f, true);
		
		tol = new Tolerance(20f,true);
		
		WindowFilter filter = new WindowFilter(6, 50);
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		PeakParameter.setPartitionNum(1);
		//PeakParameter.setGroupNum(1);
		PeakParameter.setIonCurrentRatioNum(1);
		
		SpectrumParameter.setSpecMzRangeNum(1);
		
		IonSelector ionfinder = new IonSelector(mgf, tol, filter, 8);
		
		ionfinder.train(charge, false);
		
		ArrayList<IonType> sigIons = new ArrayList<IonType>();
		for(IonType ion : ionfinder.getSigIons().keySet()){
			if(!ion.equals(IonType.NOISE)) 
				sigIons.add(ion);
		}
		
		
		FeatureSelector confinder = new FeatureSelector(mgf, sigIons, tol, pmtol).filter(filter);
		confinder.sigprob = 0.15f;
		confinder.train(charge, 0, aaSet);
		
		FeatureStatisticsTrainer contrainer = new FeatureStatisticsTrainer(mgf, ionfinder.getSigIons(), confinder.getSignificantConditions(), tol, pmtol).filter(filter);
		contrainer.train(charge, 0);
		contrainer.discardConditionsWithLessThan(100, 5000, true);
		contrainer.writeInFile(para, charge, 0, true);
		
	}
}
