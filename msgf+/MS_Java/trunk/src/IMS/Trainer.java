package IMS;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class Trainer {
	private ArrayList<Spectrum> spectra;
	private ArrayList<IonType> usedIons;
	private Tolerance tol;
	
	private Trainer(String mgf, ArrayList<IonType> usedIons, Tolerance tol) throws IOException{
		this.usedIons = usedIons;
		this.tol = tol;
		spectra = new ArrayList<Spectrum>();
		Iterator<Spectrum> iterator = new SpectraIterator(mgf, new MgfSpectrumParser());
		while(iterator.hasNext()) spectra.add(iterator.next());
		System.out.println("Total spec num: " + spectra.size());
	}
	
	static private Peptide getReversedPeptide(Peptide p){
		ArrayList<AminoAcid> rp = new ArrayList<AminoAcid>();
		for(int i=p.size()-2;i>=0;i--)
			rp.add(p.get(i));
		rp.add(p.get(p.size()-1));
		return new Peptide(rp);		
	}
	
	private HashMap<Integer, FragmentParameter> getFragmentParameters(Spectrum s, boolean isTarget){
		HashMap<Integer, FragmentParameter> paras = new HashMap<Integer, FragmentParameter>();
		Peptide peptide = s.getAnnotation();
		if(!isTarget) peptide = getReversedPeptide(peptide);
		for(int r=1;r<peptide.size();r++){
			paras.put(r, new FragmentParameter(peptide, r));
		}
		return paras;
	}
	
	private HashMap<Integer, HashMap<IonType, Float>> getIons(Spectrum s, boolean isTarget){		
		HashMap<Integer, HashMap<IonType, Float>> ions = new HashMap<Integer, HashMap<IonType, Float>>();
		Peptide peptide = s.getAnnotation();
		if(!isTarget) peptide = getReversedPeptide(peptide);
		for(int r=1;r<peptide.size();r++){
			float prefixMass = peptide.getMass(0, r);
			float suffixMass = peptide.getMass(r, peptide.size());
			HashMap<IonType, Float> fions = new HashMap<IonType, Float>();
			for(IonType ion : usedIons){
				float targetMz = 0;
				if(ion.isPrefixIon()){
					targetMz = ion.getMz(prefixMass);
				}else if(ion.isSuffixIon()){
					targetMz = ion.getMz(suffixMass);
				}else continue;
				Peak peak = s.getPeakByMass(targetMz, tol);
				if(peak != null){
					fions.put(ion, peak.getIntensity());
				}else{
					fions.put(ion, 0f);
				}
			}
			ions.put(r, fions);			
		}		
		return ions;
	}
	
	
	private void train(boolean isTarget){
		int i=0;
		for(Spectrum s : spectra){
			HashMap<Integer, FragmentParameter> paras = getFragmentParameters(s, isTarget);
			HashMap<Integer, HashMap<IonType, Float>> ions = getIons(s, isTarget);
			
			for(int r : paras.keySet()){
				IonCounter.update(paras.get(r), ions.get(r), isTarget);
				IonPairCounter.update(paras.get(r), ions.get(r), isTarget);
			}
			System.out.println("Spec num: " + i++);			
		}		
	}
	
	private void writeProbs(String file, boolean isTarget) throws IOException{
		PrintStream out = new PrintStream(file);
		out.println("##Ion Probs");
		for(IonType ion : usedIons) out.println("## " + ion);
		for(FragmentParameter para: FragmentParameter.getAllFragmentParameters()){
			out.print(para.toFileString()+"\t");
			for(IonType ion : usedIons){
				float prob = IonCounter.getProbability(para, ion, isTarget);
				out.print(prob+"\t");
			}
			out.println();
		}
		
		out.println("##Ion Pair Probs");
		
		for(FragmentParameter para: FragmentParameter.getAllFragmentParameters()){
			out.println("PARA\t" + para.toFileString());
			for(IonType ion1 : usedIons){				
				for(IonType ion2 : usedIons){
					if(ion1.equals(ion2) || ion1.getCharge()!=ion2.getCharge()) continue;
					out.println("ION1\t"+usedIons.indexOf(ion1)+"\tION2\t"+usedIons.indexOf(ion2));
					for(int r : IonPairCounter.getAllRatios()){
						float prob = IonPairCounter.getProbability(para, ion1, ion2, r, isTarget);
						out.print(r+":"+prob+"\t");
					}
					out.println();
				}
			}
		}
		
		out.close();
	}
	
	private void writeScores(String file) throws IOException{
		PrintStream out = new PrintStream(file);
		out.println("##\t1");
		for(IonType ion : usedIons) out.println("ION\t" + ion.getCharge() + "\t" + ion.getName());
		for(FragmentParameter para: FragmentParameter.getAllFragmentParameters()){
			out.print(para.toFileString()+"\t");
			for(IonType ion : usedIons){
				float probt = IonCounter.getProbability(para, ion, true);
				float probd = IonCounter.getProbability(para, ion, false);
				out.print(Math.log(probt/probd)+"\t");
			}
			out.println();
		}
		
		out.println("##\t2");
		
		for(FragmentParameter para: FragmentParameter.getAllFragmentParameters()){
			out.println("PARA\t" + para.toFileString());
			for(IonType ion1 : usedIons){				
				for(IonType ion2 : usedIons){
					if(ion1.equals(ion2) || ion1.getCharge()!=ion2.getCharge()) continue;
					out.println("ION1\t"+usedIons.indexOf(ion1)+"\tION2\t"+usedIons.indexOf(ion2));
					for(int r : IonPairCounter.getAllRatios()){
						float probt = IonPairCounter.getProbability(para, ion1, ion2, r, true);
						float probd = IonPairCounter.getProbability(para, ion1, ion2, r, false);
						out.print(r+":"+Math.log(probt/probd)+"\t");
					}
					out.println();
				}
			}
		}
		
		out.close();
	}
	
	public static void main(String[] args) throws IOException{
		String mgf = "/home/kwj/Dropbox/Training/HCD_train.mgf";
		String score = "/home/kwj/Dropbox/IMS/ScorePara.txt";
		String target = "/home/kwj/Dropbox/IMS/tragetPara.txt";
		String decoy = "/home/kwj/Dropbox/IMS/decoyPara.txt";
		
		ArrayList<IonType> usedIons = new ArrayList<IonType>();
		String[] base ={
				"y","a","b","y2","a2","b2","y3","a3","b3",
        };
        String[] extension={
        		"", "-H2O", "-NH3", 
        };
        for (int i = 0; i < base.length; i++) {
            for (int j = 0; j < extension.length; j++) {
            	usedIons.add(IonType.getIonType(base[i]+extension[j]));
            }
        }
        
        Trainer trainer = new Trainer(mgf, usedIons, new Tolerance(25f, true));
        trainer.train(true);
        trainer.writeProbs(target, true);
        trainer.train(false);        
        trainer.writeProbs(decoy, false);		
        trainer.writeScores(score);		
	}
}
