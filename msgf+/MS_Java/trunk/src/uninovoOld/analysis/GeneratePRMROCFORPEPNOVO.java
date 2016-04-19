package uninovoOld.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import msgf.Tolerance;
import msutil.Peak;
import msutil.Peptide;
import msutil.Spectrum;
import parser.BufferedLineReader;


public class GeneratePRMROCFORPEPNOVO {
	static int numNodes = 20;
	static Tolerance tol = new Tolerance(0.5f, false);
	static int charge = 2;
	
	static public void main(String[] args) throws IOException{
		String PRMfilename = "/home/kwj/Dropbox/PRM/CIDTrypsin_testPepNovoPRM.txt";
		
		if(PRMfilename.contains("HCD")) tol = new Tolerance(20f, true);
		String s;
		BufferedLineReader in = new BufferedLineReader(PRMfilename);
		HashMap<Integer, int[]> scores = new HashMap<Integer, int[]>();
		
		int[] scoreRange = new int[2];
		boolean begun = false;
		Peptide pep = null;
		int c = -1;
		Spectrum PRM = null;
		int totalt = 0, totalf = 0;
		
		while((s=in.readLine())!=null){
			if(s.isEmpty()) continue;
			//System.out.println(s);
			String[] t = s.split(" ");
			if(s.startsWith(">>")){
				
				if(begun && !PRM.isEmpty()){
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
						
					}
					PRM.setRanksOfPeaks();
				//	System.out.println(PRM);
					for(int i=0; i<PRM.size();i++){ // exclude sink source
						Peak p = PRM.get(i);
						if(p.getRank() > numNodes) continue;
						
						int score = -p.getRank();//Math.round(p.getIntensity());
						
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
				
				
				begun = true;
				
				pep = new Peptide(t[t.length-1]);
				
				continue;
			}
			if(!begun) continue;
			
			if(s.startsWith("Charge")){
				c = Integer.parseInt(t[1]);
				PRM = new Spectrum(0, c, 100);
				
				PRM.setAnnotation(pep);
				PRM.correctParentMass(pep);
				continue;
			}
			if(charge > 0 && c != charge) continue;
			
			PRM.add(new Peak(Float.parseFloat(t[0]), Float.parseFloat(t[1]), 1));
		}
		
		in.close();
		
		
		int t = 0, f =0;
		//totalt = 0;
		for(int i : scores.keySet()){
		//	totalt += scores.get(i)[0];
			totalf += scores.get(i)[1];
		}
		//System.out.println(scoreRange[0] + "\t" + scoreRange[1]);
		System.out.println("pa=[");
		for(int i=scoreRange[1]; i>=scoreRange[0]; i--){
			if(!scores.containsKey(i)) continue;
			t += scores.get(i)[0];
			f += scores.get(i)[1];
			System.out.println((float)f/totalf + "\t" + (float)t/totalt);
		}
		System.out.println("];");
		
		System.out.println("pb=[");
		float sum = 0;
		for(int i=scoreRange[1]; i>=Math.max(scoreRange[0], scoreRange[1]-50); i--){
			if(!scores.containsKey(i)) continue;
			sum += (float)scores.get(i)[0];
			System.out.println(-i + "\t" + sum/totalt);
		}
		System.out.println("];%"+totalt+"\t" +totalf);
	}
}
