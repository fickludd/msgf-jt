package uninovoOld.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class GeneratePRMROC {
	
	static int numNodes = 100;
	static Tolerance tol = new Tolerance(0.5f, false);
	static int charge = 2;
	
	static public void main(String[] args) throws IOException{
		String PRMspecfilename = "/home/kwj/Dropbox/UniNovo/PRM/HCD.mgf";
	//	if(!PRMspecfilename.contains("MSGF") && PRMspecfilename.contains("HCD"))tol = new Tolerance(20f, true);
		//System.out.println(tol);
		Iterator<Spectrum> iterator = new SpectraIterator(PRMspecfilename, new MgfSpectrumParser());
		HashMap<Integer, int[]> scores = new HashMap<Integer, int[]>();
		
		int[] scoreRange = new int[2];
		int totalt = 0, totalf = 0;
		int sn = 0;
		
		while(iterator.hasNext()){
			Spectrum PRM = iterator.next();			
			Spectrum SRM = null;
			if(PRMspecfilename.contains("MSGF")){
				SRM = iterator.next();
				for(Peak p : PRM){
					//System.out.println(PRM.getTitle()+ "\t" + p+ "\t" + SRM.getPeakByMass(PRM.getPeptideMass()-p.getMz(), tol));
					if(SRM.getPeakByMass(PRM.getPeptideMass()-p.getMz(), tol) != null)
						p.setIntensity(p.getIntensity() + SRM.getPeakByMass(PRM.getPeptideMass()-p.getMz(), tol).getIntensity());
				}
			}
			
			if(charge > 0 && PRM.getCharge() != charge) continue;
			sn++;
			
			totalt += PRM.getAnnotation().size()-1;
			
			if(PRM.get(0).getMz() == 0f){
				PRM.remove(0);
			}
			if(Math.abs(PRM.get(PRM.size()-1).getMz() - PRM.getAnnotation().getMass()) < tol.getToleranceAsDa(PRM.getAnnotation().getMass())){
				PRM.remove(PRM.size()-1);
			}
			
			
			
			ArrayList<Float> correctPRMs = new ArrayList<Float>();
			
			for(float m : PRM.getAnnotation().getPRMMasses(true, 0)){
				correctPRMs.add(m);
				//if(PRM.getPeakByMass(m, tol) == null){
				//	PRM.add(new Peak(m, -100, 1));
				//}
			}
			
			PRM.setRanksOfPeaks();
			
			for(int i=0; i<PRM.size();i++){ // exclude sink source
				Peak p = PRM.get(i);
				if(p.getRank() > numNodes) continue;
				
				int score = -p.getRank();//Math.round(p.getIntensity());
				//if(score == 29) System.out.println(PRM.getAnnotationStr()+"\t" + p.getIntensity());
				boolean isTrue = false;
				for(float cm : correctPRMs){
					if(Math.abs(cm - p.getMz()) <= tol.getToleranceAsDa(Math.max(cm, PRM.getAnnotation().getMass()-cm))){
						isTrue = true;
						break;
					}
				}
				
				if(!scores.containsKey(score)){
					scores.put(score, new int[2]);
				}
				
				if(isTrue){
					scores.get(score)[0] ++;
				}else{
					scores.get(score)[1] ++;
				}
				scoreRange[0] = Math.min(scoreRange[0], score);
				scoreRange[1] = Math.max(scoreRange[1], score);
			}
			
			
		}
		
	
		int t = 0, f =0;
		int ttotalt=0;
		for(int i : scores.keySet()){
			ttotalt += scores.get(i)[0];
			totalf += scores.get(i)[1];
		}
		//System.out.println(scoreRange[0] + "\t" + scoreRange[1]);
		String prefix = "u";
		if(PRMspecfilename.contains("MSGF")) prefix = "m";
		System.out.println(prefix+"a=[");
		for(int i=scoreRange[1]; i>=scoreRange[0]; i--){
			if(!scores.containsKey(i)) continue;
			t += scores.get(i)[0];
			f += scores.get(i)[1];
			System.out.println((float)f/totalf + "\t" + (float)t/ttotalt);
		}
		System.out.println("];");
		
		float sum = 0;
		System.out.println(prefix+"b=[");
		for(int i=scoreRange[1]; i>=Math.max(scoreRange[0], scoreRange[1]-50); i--){
			if(!scores.containsKey(i)) continue;
			sum += (float)scores.get(i)[0];
			System.out.println(-i + "\t" + sum/totalt);///(scores.get(i)[0] + scores.get(i)[1])
		}
		System.out.println("];%"+totalt+"\t" +totalf + "\t" + sn);
		
	}
}
