package uninovoOld.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import msutil.Enzyme;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;

public class MakeTrainingAndTestMgf {
	static int testSetSizePerCharge = 1000;
	static int maxNumSpecFortheSamePeptide = 3;
	
	
	static public void main(String[] args){
		String inputmgf = "/home/kwj/workspace/inputs/MultiEnzymes/ETD_IT_AspN/";//"/home/kwj/workspace/inputs/HeckRevision/AnnotatedSpectra/CID_LysN_Confident.mgf";
		String inputmgf2 = null;// "/home/kwj/workspace/inputs/CIDETDPairs/annotated/ETDTrypsin_PAIRED.mgf";
		Enzyme enzyme = Enzyme.AspN;
		
		try {
			
			if(new File(inputmgf).isDirectory()){
				String newfile = inputmgf+"/merged.mgf";
				PrintStream ps = new PrintStream(newfile);
				
				for(File f : new File(inputmgf).listFiles()){
					if(f.getAbsolutePath().endsWith(".mgf") && !f.getAbsolutePath().endsWith("merged.mgf")){
						BufferedLineReader in = new BufferedLineReader(f.getAbsolutePath());
						System.out.println(f);
						String s;
						while((s = in.readLine()) != null){
							ps.println(s);
						}
						in.close();
					}
				}
				ps.close();
				inputmgf = newfile;
			}
			
			Iterator<Spectrum> iterator = new SpectraIterator(inputmgf, new MgfSpectrumParser());
			Iterator<Spectrum> iterator2 = null;
			if(inputmgf2 != null) iterator2 = new SpectraIterator(inputmgf2, new MgfSpectrumParser());
			HashMap<String, ArrayList<Spectrum>> specs = new HashMap<String, ArrayList<Spectrum>> ();
			HashMap<String, ArrayList<Spectrum>> specs2 = null;
			if(inputmgf2 != null) specs2 = new HashMap<String, ArrayList<Spectrum>> ();
			PrintStream test = new PrintStream(inputmgf.substring(0, inputmgf.lastIndexOf(".")) + "_test.mgf");
			PrintStream test2 = null;
			if(inputmgf2!=null) test2 = new PrintStream(inputmgf2.substring(0, inputmgf.lastIndexOf(".")) + "_test.mgf");
			PrintStream train = new PrintStream(inputmgf.substring(0, inputmgf.lastIndexOf(".")) + "_train.mgf");
			
			
			int[] chargeRange = new int[2];
			chargeRange[0] = Integer.MAX_VALUE;
			
			Random r = new Random();
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				Spectrum spec2 = null;
				if(iterator2 != null){
					spec2 = iterator2.next();
					//System.out.println(spec.getAnnotationStr() + "\t" + spec2.getAnnotationStr());
					assert(spec.getParentMass() == spec2.getParentMass());
				}
				
				if(spec.getAnnotation() == null) continue;
				if(enzyme!= null && !enzyme.isCleaved(spec.getAnnotation())) continue;
				if(spec.getAnnotation().isModified()) continue;
				if(spec.getAnnotationStr().isEmpty()) continue;
				if(Math.abs(spec.getAnnotation().getMass() - spec.getPeptideMass()) > 1) continue;
				
				String pep = spec.getAnnotationStr();
				pep = spec.getCharge() + pep;
				
				chargeRange[0] = Math.min(chargeRange[0], spec.getCharge());
				chargeRange[1] = Math.max(chargeRange[1], spec.getCharge());
				
				int n = 0;
				
				if(!specs.containsKey(pep)){
					specs.put(pep, new ArrayList<Spectrum>());
					if(specs2 != null) specs2.put(pep, new ArrayList<Spectrum>());
				}
				
				ArrayList<Spectrum> ss = specs.get(pep);
				ArrayList<Spectrum> ss2 = null;
				if(specs2 != null) ss2 = specs2.get(pep);
				
				if(ss.size() < maxNumSpecFortheSamePeptide){
					ss.add(spec);
					if(ss2 != null) ss2.add(spec2);
				}else{
					n = r.nextInt(maxNumSpecFortheSamePeptide+1);
					if(n < maxNumSpecFortheSamePeptide){
						ss.remove(n);
						ss.add(spec);
						if(ss2 != null){
							ss2.remove(n);
							ss2.add(spec2);
						}
					}
				}				
			}
			
			System.out.println(specs.size());
			ArrayList<Spectrum> testSet = new ArrayList<Spectrum>();
			HashSet<String> testPeptides = new HashSet<String>();
			
			int counter = 0;
			for(int charge = chargeRange[0]; charge <= chargeRange[1]; ){
				int counter2 = 0;
				for(String pep : specs.keySet()){
					if(pep.startsWith(new Integer(charge).toString())){
						counter2 ++;	
						int n = r.nextInt(specs.size()/testSetSizePerCharge);
						if(specs.get(pep).isEmpty()) continue;
						Spectrum s = specs.get(pep).get(0);
						Spectrum s2 = null;
						if(specs2 != null) s2 = specs2.get(pep).get(0);
						
						
						if(n==0 && !testSet.contains(s)){
							testSet.add(s);
							
							testPeptides.add(pep);
							s.outputMgf(test,true);
							if(s2 != null){
								s2.outputMgf(test2, true);
								System.out.println(s.getParentMass() + "\t" + s2.getParentMass());
								assert(s.getParentMass() == s2.getParentMass());
							}
							specs.get(pep).remove(0);
							if(specs2 != null)specs2.get(pep).remove(0);
							counter++;
							System.out.println(charge + "\t" + counter);
							if(counter >= testSetSizePerCharge){
								counter = 0;
								charge++;
								break;
							}
						}
					}
				}	
				if(counter2 < testSetSizePerCharge){
					if(counter > 0){
						charge ++;
						counter = 0;
					}
					
				}
			}
			System.out.println("Test set done");
			int cn=0;
			for(String pep : specs.keySet()){
				if(testPeptides.contains(pep)) continue;
				ArrayList<Spectrum> ss = specs.get(pep);
				for(Spectrum s : ss){
					s.outputMgf(train);
					if(s.getCharge() == 3) cn ++;
				}
			}
			
			System.out.println("Training set done" +cn);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}
