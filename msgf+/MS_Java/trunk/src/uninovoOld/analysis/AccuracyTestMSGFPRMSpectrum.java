package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;
import uninovoOld.DenovoReconstruction;
import uninovoOld.SpectrumGraph;
import uninovoOld.UniNovo;

public class AccuracyTestMSGFPRMSpectrum {
	static String output="";
	
	static public void main(String[] args) throws IOException{
		String prmSpectrumFile = "/home/kwj/Dropbox/Test/ETDAspN_testPRMMSGF.mgf";
		int charge = 4;
		Enzyme enzyme = Enzyme.AspN;//Enzyme.TRYPSIN;
		
		int numPerLength = 20;
		
		int minNumPerLength = numPerLength;
		int minLength = 5;
		int maxLength = 30;
		int maxGapNum = 2;
		int maxPTMNum = 0;
		boolean accuracyFirst = false;
		float[] gapnum = new float[50];
		int[] numcorrect = new int[50];
		float[] avglength = new float[50];
		float totalcorrectlength = 0;
		int totalnumcorrect = 0;
		int sn = 0;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		Tolerance tol = new Tolerance(0.5f, false);
		if(prmSpectrumFile.contains("HCD") && !prmSpectrumFile.contains("MSGF")) tol = new Tolerance(20f, true);
		Tolerance pmtol = new Tolerance(20f, false);
		if(prmSpectrumFile.contains("MSGF")) pmtol = new Tolerance(0.5f, false);
		Iterator<Spectrum> iterator = new SpectraIterator(prmSpectrumFile, new MgfSpectrumParser());
		
		if(output != null && prmSpectrumFile.contains("MSGF")){
			output = prmSpectrumFile + ".out"+numPerLength+"_"+ charge;
		}
		
		PrintStream outfile = null;
		if(output != null) outfile = new PrintStream(output);
		
		
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			Spectrum SRM = null;
			if(prmSpectrumFile.contains("MSGF")){				
				SRM = iterator.next();
				for(Peak p : s){
					if(SRM.getPeakByMass(s.getPeptideMass()-p.getMz(), tol) != null)
						p.setIntensity(p.getIntensity() + SRM.getPeakByMass(s.getPeptideMass()-p.getMz(), tol).getIntensity());
				}
			}
			
			if(enzyme!= null && !enzyme.isCleaved(s.getAnnotation())) continue;
			
			//if(s.getScanNum() < 1000) continue;
			if(charge > 0 && s.getCharge() != charge) continue;
			if(s.getAnnotation().isModified()) continue;
			
			sn++;
			SpectrumGraph g = SpectrumGraph.getSpectrumGraphFromPRMSpectrum(s, tol, pmtol, aaSet, null);
			
			UniNovo mg = new UniNovo(g);
			//System.out.println(g.getNodes());
			ArrayList<DenovoReconstruction> out = new ArrayList<DenovoReconstruction>();
			mg.getOutput(out, numPerLength, minNumPerLength, minLength, maxLength, 0, maxGapNum, maxPTMNum, accuracyFirst);
		
			int size = 0;
			float correctlen = 0;
			boolean isCorrect = false;
			System.out.println(s.getScanNum() + "\t" + s.getAnnotationStr() + "\t" + s.getCharge());
			
			for(DenovoReconstruction gp : out){
				System.out.print(gp + " : " + Math.round(gp.getAccuracy()*100)+"% ");
						
				if(gp.isCorrect(s.getAnnotation())){
					isCorrect = true;
					System.out.println(" true*\t" + gp.length()  + "\t" + gp.getLR());
					if(correctlen == 0) correctlen = gp.length();
					//			lengthLRaccuracy[gp.length()][Math.max(0, Math.round(gp.getLR()))]++;
				}else{
					System.out.println(" false\t" + gp.length() + "\t" +  gp.getLR());						
				}
	
				size++;
				gapnum[0] += gp.getGapNum();
		
			}
			
			if(isCorrect){
				avglength[0] += correctlen;
				totalcorrectlength += correctlen;
				numcorrect[0]++;
				totalnumcorrect++;
				if(output != null)  outfile.println(s.getCharge()+"\t"+s.getAnnotationStr() + "\t" + correctlen);

			}
		}
		
		if(output != null)  outfile.println("###\t" + 0 + "\t" + 0 + "\t" + sn + "\t" + totalnumcorrect);
		if(output != null)  outfile.close();
		
		System.out.println("accuracy: " + ((float)totalnumcorrect/sn));
	//	System.out.println("total accuracy: " + ((float)totalnumcorrect/sn));
		System.out.println("length: " + totalcorrectlength/totalnumcorrect);
		//System.out.println("correct cut number: " + correctCutLen + "\t" + incorrectCutLen);
	}
}
