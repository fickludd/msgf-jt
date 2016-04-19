package uninovoOld;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import uninovoOld.features.Feature;
import uninovoOld.features.IntensityFeature;
import uninovoOld.features.LinkingFeature;
import uninovoOld.features.OffsetFeature;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;

public class IPDGenerator {
	
	static private HashMap<Integer, HashMap<Integer, HashMap<Integer, HashMap<IonType, Float>>>> ionChargeMap = null; // type , iteration, ion, intensity
	static private HashMap<Integer, HashMap<Integer, HashMap<Integer, ArrayList<IonType>>>> sortedIonChargeMap = null;
	static  private HashMap<Integer, HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>> independentConditionMap = null;
	static  private HashMap<Integer, HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>>>> rankIonProbMap = null;
	
	static private Tolerance[] tols;
	static private Tolerance[] pmtols;
	static private boolean isForPMCorrection = false;
	
	private String parafile;
	
	private AminoAcidSet aaSet;
	//private float minKLDivergence = 0;//1e-2f;
	private WindowFilter filter = null;
	private int typeIndex = 0;
	private int maxIterationNum = 0;
	private int maxCharge = 100;
	private int minCharge = 2;
	
	static public void clear(){
		ionChargeMap = null;
		sortedIonChargeMap = null; 
		independentConditionMap = null;
		rankIonProbMap = null;
	}
	
	
	public IPDGenerator(String parafile, AminoAcidSet aaSet, int typeIndex){
		this.parafile = parafile;
		this.aaSet = aaSet;
		this.typeIndex = typeIndex;
		readFromFile();
		maxCharge = getChargeRange().get(1);
		minCharge = getChargeRange().get(0);
	}
	
	static public int getMaxRank(Spectrum spec){ return Math.round(spec.getPeptideMass()/121.6f * 20);}//Integer.MAX_VALUE;}
	
	public IPDGenerator filter(WindowFilter b) {filter = b; return this;}
	
	public AminoAcidSet getAASet() {return aaSet;}
	
	public String getParaFile() {return parafile;}
	
	public ArrayList<Integer> getChargeRange(){
		ArrayList<Integer> r = new ArrayList<Integer>();
		int min = 1000, max = 0;
		for(int c : ionChargeMap.get(typeIndex).get(0).keySet()){
			min = Math.min(min, c);
			max = Math.max(max, c);
		}
		r.add(min); r.add(max);
		return r;
	}
	
	public ArrayList<IonType> getSigIonsOrderedByIntensityWithOutNoiseIon(int charge) {
		HashMap<Integer, HashMap<Integer, HashMap<IonType, Float>>> ionMap = ionChargeMap.get(typeIndex);
		if(ionMap == null) return new ArrayList<IonType>();
		
		int max = 0;
		for(int n : ionMap.keySet()) max = Math.max(max, n);
		
		return getSigIonsOrderedByIntensityWithOutNoiseIon(charge, max);	
	}
	
	
	
	public ArrayList<IonType> getSigIonsOrderedByIntensityWithOutNoiseIon(int charge, int iterationNumber) {
		charge = Math.min(charge, maxCharge);
		charge = Math.max(charge, minCharge);
		
		if(ionChargeMap.get(typeIndex) == null || ionChargeMap.get(typeIndex).get(iterationNumber) == null) return new ArrayList<IonType>();
		
		if(sortedIonChargeMap == null) sortedIonChargeMap = new HashMap<Integer, HashMap<Integer, HashMap<Integer, ArrayList<IonType>>>>();
		
		if(!sortedIonChargeMap.containsKey(typeIndex)){
			sortedIonChargeMap.put(typeIndex, new HashMap<Integer, HashMap<Integer, ArrayList<IonType>>>());
		}
		
		HashMap<Integer, HashMap<Integer, ArrayList<IonType>>>
			subsortedIonChargeMap = sortedIonChargeMap.get(typeIndex);
		
		if(!subsortedIonChargeMap.containsKey(iterationNumber)){
			subsortedIonChargeMap.put(iterationNumber, new HashMap<Integer, ArrayList<IonType>>());
		}
		
		HashMap<Integer, ArrayList<IonType>>
			subsubsortedIonChargeMap = subsortedIonChargeMap.get(iterationNumber);
		
		if(!subsubsortedIonChargeMap.containsKey(charge))
			subsubsortedIonChargeMap.put(charge, new ArrayList<IonType>());
		
		
		ArrayList<IonType> ret = subsubsortedIonChargeMap.get(charge);
		
		if(ret != null && !ret.isEmpty()) return ret;
		
		ArrayList<Float> intensities = new ArrayList<Float>();
		//
		//System.out.println(ionChargeMap.get(typeIndex).get(iterationNumber));
		for(IonType ion : ionChargeMap.get(typeIndex).get(iterationNumber).get(charge).keySet()){
			intensities.add(ionChargeMap.get(typeIndex).get(iterationNumber).get(charge).get(ion));
		}
		
		Collections.sort(intensities);
		
		for(int i = intensities.size()-1; i>=0; i--){
			for(IonType ion : ionChargeMap.get(typeIndex).get(iterationNumber).get(charge).keySet()){
				if(intensities.get(i) == ionChargeMap.get(typeIndex).get(iterationNumber).get(charge).get(ion)){
					ret.add(ion);
				}
			}
		//	if(ret.size() > num) break;
		}
		ret.remove(IonType.NOISE);
		//System.out.println(charge + "\t" + iterationNumber + "\t" + ret);
		return ret;
	}
	
	
	private void updateSpectrum(Spectrum spec, HashMap<Peak, HashMap<IonType, Float>> profile, int iterationNumber, boolean isLastIteration, int emphasizeIthIon){
		 float high = 0, low = 0;
		 //int charge = Math.min(maxCharge, spec.getCharge());
		// charge = Math.max(minCharge, charge);
		// spec.setRanksOfPeaks();
		 
		//ArrayList<IonType> sigIons = new ArrayList<IonType>();
		
		//sigIons.addAll(getSigIonsOrderedByIntensityWithOutNoiseIon(charge, iterationNumber));
		 SpectrumParameter spar = new SpectrumParameter(spec);
		 
		//System.out.println(sigIons);
	//	if(!sigIons.contains(IonType.NOISE))sigIons.add(IonType.NOISE);
		HashMap<IonType, Float> ionfactorMap = ionChargeMap.get(typeIndex).get(iterationNumber).get(spar.getSpecCharge());
			
		 for(Peak p : spec){
			 if(p.getRank() == 1){
				 high = p.getIntensity();
				 break;
			 }
		}
		
	//	 if(tol.getToleranceAsDa(500)*2 >= 0.5f)
		 low = high*1e-2f; // bad for high precision, good for low precision.
		 
		 for(Peak p : profile.keySet()){
			HashMap<IonType, Float> probs = profile.get(p);
			
			//int pp = new PeakParameter(p, spec, iterationNumber, minCharge, maxCharge).getBasePeakPartitionNum();
			
			float intensity = 0; 
				
			float factor = 1;
			
			ArrayList<IonType> ions = getSigIonsOrderedByIntensityWithOutNoiseIon(spec.getCharge(), iterationNumber);
			
			for(int i=0; i<ions.size(); i++){
				IonType ion = ions.get(i);
				if(ion instanceof IonType.PrecursorIon) continue;
				if(ion.equals(IonType.NOISE)) continue;
				
				if(isLastIteration && emphasizeIthIon >=0 && i != emphasizeIthIon) continue; //
				
				//if(isLastIteration)
					factor = ionfactorMap.get(ion);
				//else{
				//	factor = (float) (1 -  Math.pow((float)i/sigIons.size(),2));
				//}
				//	if(isLastIteration) factor = 10000000;
				intensity += probs.get(ion) *  factor;

			}
			intensity = intensity * (high - low) + low;
			p.setIntensity(intensity);
		 }
		 
		 /*for(int i=0; i< spec.size(); i++){
			 Peak p = spec.get(i);
			 if(p.getIntensity() == 0){
				 System.out.println("******");
				 spec.remove(i--);
			 }
		}*/
	}
	
	private HashMap<Peak, HashMap<IonType, Float>> getIPDForEachIteration(Spectrum spec, int maxRank, int iterationNum){
		HashMap<Peak, HashMap<IonType, Float>> profile = new HashMap<Peak, HashMap<IonType, Float>>();	
		spec.setRanksOfPeaks();
		
		Spectrum filteredspec = spec;
		//int charge = Math.min(maxCharge, spec.getCharge());
	//	charge = Math.max(minCharge, charge);
		
		if(filter != null) filteredspec = filter.apply(spec);
						
		SpectrumParameter spar = new SpectrumParameter(spec);
		
		HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> sigconMap = independentConditionMap.get(typeIndex).get(iterationNum).get(spar);
		//HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>> srankIonMap = rankIonProbMap.get(typeIndex).get(iterationNum).get(spar);
	
		if(sigconMap == null || sigconMap.isEmpty()) return profile;
		
		PeakParameter.refresh();
		
	//	if(iterationNum > 0 ) maxRank = (int)(maxRank * 0.5);
		for(Peak bp : spec){
			if(bp.getRank() > maxRank) continue;
			PeakParameter ppar = new PeakParameter(bp, spec, iterationNum);
			
			if(ppar.getBasePeakGroupNum() > PeakParameter.getMaxGroupNum()) continue;

			HashMap<Integer, ArrayList<Feature>> consMap = sigconMap.get(ppar);
			
			/*HashMap<Integer, HashMap<IonType, Float>> ssrankIonMap = srankIonMap.get(ppar);
			HashMap<IonType, Float> sssrankIonMap = null;
			if(ssrankIonMap != null){
				sssrankIonMap = ssrankIonMap.get(bp.getRank());
				if(sssrankIonMap == null){
					int mind = Integer.MAX_VALUE; 
				
					for(int r : ssrankIonMap.keySet()){
						int d = Math.abs(bp.getRank() - r);
						if(mind > d){
							mind = d;
							sssrankIonMap = ssrankIonMap.get(r);
						}
					}
				}
			}*/
		//	sssrankIonMap = null;
			ArrayList<Feature> cons = new ArrayList<Feature>();
			
			for(int setnum : consMap.keySet()){	
				ArrayList<Feature> features = consMap.get(setnum);
				
				assert(!features.isEmpty());
				
				for(int i = 0; i< features.size(); i++){
					Feature con = features.get(i);
					if(isForPMCorrection()){						
						if(con instanceof OffsetFeature){
							OffsetFeature ocon = (OffsetFeature)con;
							if(!ocon.isComplementary()) continue;
						}else if(!(con instanceof IntensityFeature))
							continue;
					}
					if(con.holdsFor(bp, filteredspec, spar, ppar, tols[typeIndex], pmtols[typeIndex], iterationNum)){
						cons.add(con);
						break;
					}			//if(i>50)break;
				}
			}
			
			if(!cons.isEmpty())
				profile.put(bp, Feature.getPeakIonProbabilityDistribution(bp, cons, getSigIonsOrderedByIntensityWithOutNoiseIon(spar.getSpecCharge(), iterationNum), null, spec));
			
		}
		return profile;
	}
	
	private void updateIndependentConditionMap(Feature con, int setNum){
		int iterationNum = con.getIterationNum();
		SpectrumParameter spar = con.getSpectrumParameter();
		PeakParameter ppar = con.getBasePeakParameter();
		
		if(!independentConditionMap.containsKey(typeIndex))
			independentConditionMap.put(typeIndex, new HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>());
			
		if(!independentConditionMap.get(typeIndex).containsKey(iterationNum))
			independentConditionMap.get(typeIndex).put(iterationNum, new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>());
		HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>
			s1 = independentConditionMap.get(typeIndex).get(iterationNum);
		
		if(!s1.containsKey(spar))
			s1.put(spar, new HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>());
		HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>
			s2 = s1.get(spar);
	
		if(!s2.containsKey(ppar))
			s2.put(ppar, new HashMap<Integer, ArrayList<Feature>>());
		HashMap<Integer, ArrayList<Feature>> s3 = s2.get(ppar);
	
		if(!s3.containsKey(setNum))
			s3.put(setNum, new ArrayList<Feature>());
		ArrayList<Feature> cons = s3.get(setNum);
		//
		int index = Collections.binarySearch(cons, con, Collections.reverseOrder());
		if(index < 0) index = -index-1;
		//if(Feature.getKLDivergenceFromNullCondition(ion, con) > minKLDivergence) 
		cons.add(index, con);
	//	Collections.sort(cons, new ConditionComparator(ion));
	
	}
	
	public int readFromFile(){
		if(getMaxIterationNum() > 0) return getMaxIterationNum();

		BufferedLineReader in;
		int maxIterationNum = 0;
		System.out.println("Reading parameters for type index " + typeIndex + " : " +  parafile);
	//	System.out.println("Reading parameter file... " + parafile);
		try {
			in = new BufferedLineReader(parafile);
			String s;
			int mode = -1;
			int setNum = -1;
			int chargeForIon =  -1;
			int iterationNumberForIon = -1;
		//	SpectrumParameter spar = null;

			HashMap<IonType, Float> ionProbMap = null;
			Feature con = null;
			Feature nullCon = null;
			SpectrumParameter spar = null;
			PeakParameter ppar = null;
			
			HashMap<Integer, HashMap<Integer, ArrayList<IonType>>>
				sigIonMap = new HashMap<Integer, HashMap<Integer, ArrayList<IonType>>>();
			
			while((s = in.readLine()) != null){
			
				if(s.startsWith("#ION\t")){
					chargeForIon = Integer.parseInt(s.split("\t")[1]);
					iterationNumberForIon = Integer.parseInt(s.split("\t")[2]);
					if(ionChargeMap == null)
						ionChargeMap = new HashMap<Integer, HashMap<Integer, HashMap<Integer, HashMap<IonType, Float>>>>();
					if(!ionChargeMap.containsKey(typeIndex))
						ionChargeMap.put(typeIndex, new HashMap<Integer, HashMap<Integer, HashMap<IonType, Float>>>());

					mode = 0; continue;
				}
				if(s.equals("#FEATURES")){
					if(independentConditionMap == null)
						independentConditionMap = new HashMap<Integer, HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>>();
					if(!independentConditionMap.containsKey(typeIndex))
						independentConditionMap.put(typeIndex, new HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>());

					mode = 1; 
					continue;
				}
				if(s.startsWith("#SPECPARTITION\t")){
					chargeForIon = Integer.parseInt(s.split("\t")[1]);
					mode = 2;
					continue;
				}
				if(s.startsWith("#RANKPROB")){
					iterationNumberForIon = Integer.parseInt(s.split("\t")[1]);
					chargeForIon = Integer.parseInt(s.split("\t")[2]);
					if(rankIonProbMap == null)
						rankIonProbMap = new HashMap<Integer, HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>>>>();
					if(!rankIonProbMap.containsKey(typeIndex))
						rankIonProbMap.put(typeIndex, new HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>>>());
					if(!rankIonProbMap.get(typeIndex).containsKey(iterationNumberForIon))
						rankIonProbMap.get(typeIndex).put(iterationNumberForIon, new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>>());
					
					mode = 3;
					continue;
				}
				
				if(s.startsWith("#PMTOL")){
					if(pmtols == null) pmtols = new Tolerance[100];
					pmtols[typeIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
					mode = -1;
					continue;
				}
				if(s.startsWith("#TOL")){
					if(tols == null) tols = new Tolerance[100];
					tols[typeIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
					mode = -1;
					continue;
				}
				if(s.startsWith("#MAXSPECMZ")){
					SpectrumParameter.setMaxSpecMzRange(Integer.parseInt(s.split("\t")[1]), Integer.parseInt(s.split("\t")[2]));
					mode = -1;
					continue;
				}
				if(s.startsWith("#############################################################")){
					mode = Integer.MAX_VALUE; continue;
				}
				if(s.startsWith("##NODE&EDGEPARAMETERS##")){
					mode = Integer.MAX_VALUE; 
					continue;
				}
				
				
				if(mode == 0){
					String[] token = s.split("\t");
					
					if(s.startsWith("#PARS")){
						//spar = SpectrumParameter.parseFileString(token[1]);
					//	pp = Integer.parseInt(token[1]);
						continue;
					}
					
					
					
					if(!ionChargeMap.get(typeIndex).containsKey(iterationNumberForIon))
						ionChargeMap.get(typeIndex).put(iterationNumberForIon, new HashMap<Integer, HashMap<IonType, Float>>());
				
					HashMap<Integer, HashMap<IonType, Float>> ionMap = ionChargeMap.get(typeIndex).get(iterationNumberForIon);
					
					if(!ionMap.containsKey(chargeForIon)){
						ionMap.put(chargeForIon, new HashMap<IonType, Float>());
					}
					
					HashMap<IonType, Float> ions = ionMap.get(chargeForIon);
				
					if(!token[0].startsWith("n"))
						ions.put(IonType.getIonType(token[0]), Float.parseFloat(token[1]));
					else ions.put(IonType.NOISE, Float.parseFloat(token[1]));
								
					if(!sigIonMap.containsKey(iterationNumberForIon))
						sigIonMap.put(iterationNumberForIon, new HashMap<Integer, ArrayList<IonType>>());
					
					HashMap<Integer, ArrayList<IonType>> imap = sigIonMap.get(iterationNumberForIon);
					
					if(!imap.containsKey(chargeForIon))
						imap.put(chargeForIon, new ArrayList<IonType>());
					
					ArrayList<IonType> is = imap.get(chargeForIon);
					
					is.add(IonType.getIonType(token[0]));

				}else if(mode == 1){
					if(s.startsWith("#")){
						setNum = Integer.parseInt(s.substring(1));
					}else if(s.startsWith("N")){
						con = IntensityFeature.parseFileString(s);
						nullCon = con;
					}else if(s.startsWith("G")){
						con = OffsetFeature.parseFileString(s);
					}else if(s.startsWith("B")){
						con = LinkingFeature.parseFileString(s, aaSet);
				//	}else if(s.startsWith("D")){
					//	con = DensityFeature.parseFileString(s);
					}else{
						ionProbMap = new HashMap<IonType, Float>();
						String[] token = s.split("\t");
						
						ArrayList<IonType> ions = sigIonMap.get(iterationNumberForIon).get(con.getSpectrumParameter().getSpecCharge());
						for(int i = 0; i < Math.min(token.length, ions.size()); i++){
							String t = token[i];
							IonType ion = ions.get(i);
							ionProbMap.put(ion, Float.parseFloat(t));
						}
						
						if(ions.size() < token.length)
							ionProbMap.put(IonType.NOISE, Float.parseFloat(token[token.length-1]));
						else{
							float sum = 1;
							for(IonType ion : ionProbMap.keySet()){
								sum -= ionProbMap.get(ion); 
							}
							sum = Math.max(sum, 0);
							ionProbMap.put(IonType.NOISE, sum);
						}
						
						con.setIonProbMap(ionProbMap);
						con.registerNullCondition(nullCon);
						updateIndependentConditionMap(con, setNum);
						maxIterationNum = Math.max(con.getIterationNum() + 1, maxIterationNum);
						
					}
				}else if(mode == 2){
					if(s.equals("null"))continue;
					String[] token = s.split("\t");
					float[] mzs = new float[token.length];
					for(int i=0; i<mzs.length;i++){
						mzs[i] = Float.parseFloat(token[i]);
					}
					SpectrumParameter.setPartitionMzs(mzs, chargeForIon);
					//
				}else if(mode == 3){
					HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>> sr = rankIonProbMap.get(typeIndex).get(iterationNumberForIon);
					if(s.startsWith("#SPAR")){
						spar = SpectrumParameter.parseSpectrumParameter(s.split("\t")[1]);
						if(!sr.containsKey(spar))
							sr.put(spar, new HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>());
						continue;
					}
					if(s.startsWith("#PPAR")){
						ppar = PeakParameter.parsePeakParameter(s.split("\t")[1]);
						if(!sr.get(spar).containsKey(ppar))
							sr.get(spar).put(ppar, new  HashMap<Integer, HashMap<IonType, Float>>());
						continue;
					}
					HashMap<Integer, HashMap<IonType, Float>> ssr = sr.get(spar).get(ppar);
					String[] token = s.split("\t");
					int k = Integer.parseInt(token[0]);
					if(!ssr.containsKey(k))
						ssr.put(k, new HashMap<IonType, Float>());
					HashMap<IonType, Float> sssr = ssr.get(k);
					
					for(int i=1; i<token.length; i++){
						float p = Float.parseFloat(token[i]);
						IonType ion = sigIonMap.get(iterationNumberForIon).get(chargeForIon).get(i-1);
						sssr.put(ion, p);
					}
				}
			}
			/*for(int k : ionChargeMap.keySet()){
				System.out.println(ionChargeMap.get(k));
				System.out.println(sigIonMap.get(k));
			}*/
			in.close();
			//for(int c : ionChargeMap.keySet())
			//	System.out.println(ionChargeMap.get(c));
			
			SpectrumGraphComponent.read(parafile, typeIndex);
		}catch (IOException e) {
			System.exit(1);
			e.printStackTrace();
		}
		setMaxIterationNum(maxIterationNum);
		
		System.out.println("Iteration Number : " + maxIterationNum);
		return maxIterationNum;
	}
	
	// Spectrum spec is updated
	HashMap<Peak, HashMap<IonType, Float>> getIPD(Spectrum s, int maxIterationNum){ 
		HashMap<Peak, HashMap<IonType, Float>> profile = null;
		 
		maxIterationNum = Math.min(this.maxIterationNum, readFromFile());
		if(isForPMCorrection()) maxIterationNum = 1;
		
		int maxRank = getMaxRank(s);//250;
		
		/*spec.setRanksOfPeaks();
		Spectrum s = spec.getCloneWithoutPeakList();
		int maxRank = getMaxRank(s);//250;
		//int maxRank = Math.round(spec.getPeptideMass()/1000*150);
		
		for(Peak p : spec){ 
			if(p.getRank() < maxRank) // TODO
				s.add(p.clone());//TODO erase in near future...??
		}
		*/
		for(int iterationNum = 0; iterationNum < maxIterationNum; iterationNum++){
			 profile = getIPDForEachIteration(s, maxRank, iterationNum);
			 updateSpectrum(s, profile, iterationNum, iterationNum == maxIterationNum-1, -1);
		}
		
		return profile;	
	}
	
	public void run(String specfilename, String outfilename, int specCharge,  int maxIterationNum, int iterationStartNum, boolean train, int emphasizeIthIon){
		
		long time = System.currentTimeMillis();
		
		if(train) setMaxIterationNum(0);
		maxIterationNum = Math.min(maxIterationNum, readFromFile());
		//System.out.println(maxIterationNum);
		if(train) setMaxIterationNum(maxIterationNum);
		
		int sn = 0;
		try {
			PrintStream out = null;
			//PrintStream out2 = null;
			if(outfilename != null){
				out = new PrintStream(outfilename);
				//out2 = new PrintStream(outfilename+"_10.mgf");
			}
			PrintStream matlabout = null;
			
			if(!train) matlabout =  new PrintStream(outfilename+"MATLAB");
			
		//	matlabout = null;
			
			Iterator<Spectrum> iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0;
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation() != null && spec.getAnnotation().isModified()) continue;
				Spectrum originalSpec = spec.getCloneWithoutPeakList();
				
				for(Peak p : spec){
					originalSpec.add(p.clone());
				}
				
				int origCharge = spec.getCharge();
				if(specCharge == 0){
					int charge = Math.max(spec.getCharge(), this.getMinCharge());
					
					charge = Math.min(charge, this.getMaxCharge());
					
					if(charge != spec.getCharge()){
	
						float parentMass= spec.getParentMass();
						spec.setCharge(charge);
						spec.correctParentMass(parentMass);
	
					}
				}
				//spec.correctParentMass(spec.getParentMass()-1);//TODO
				sn++;
				//if(sn >  1000)break;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("running " +": " + sn);
				}
				
				if(matlabout != null){
					matlabout.println("%scanID: "+ spec.getScanNum());
					matlabout.println("PEPTIDE='"+spec.getAnnotationStr()+"';");
					matlabout.println("PM="+ spec.getPeptideMass()+";");
					matlabout.print("%IONS=[");
					for(IonType ion : this.getSigIonsOrderedByIntensityWithOutNoiseIon(spec.getCharge(), maxIterationNum)){
						matlabout.println("%'" + ion+"'");
					}
					matlabout.println("%];");
					
					
					matlabout.println("ORIGSPEC=[");
					for(Peak p : spec){
						matlabout.println(p.getMz() + "\t" + p.getIntensity());
					}
					matlabout.println("];");
				}
				
				spec.setRanksOfPeaks();
				int maxRank = Math.round(spec.getPeptideMass()/1000*150);
				
				for(int l=0; l< spec.size();l++){ 
					Peak p = spec.get(l);
					if(p.getRank() >= maxRank){ // TODO
						spec.remove(p);
						l--;
					}
				}
				HashMap<Peak, HashMap<IonType, Float>> profile = null;
				for(int iterationNum = iterationStartNum; iterationNum < maxIterationNum; iterationNum++){
					 profile = getIPDForEachIteration(spec, maxRank, iterationNum);
					 updateSpectrum(spec, profile, iterationNum, (!train) && iterationNum == maxIterationNum-1, emphasizeIthIon);
				}
				
				if(outfilename != null){
					if(origCharge != spec.getCharge()){
					
						float parentMass= spec.getParentMass();
						spec.setCharge(origCharge);
						spec.correctParentMass(parentMass);
		
					}
					
					spec.setRanksOfPeaks();
					for(Peak p : spec){
						if(p.getRank()<=30){
							int r = p.getRank();
							float mz = p.getMz();
							
							originalSpec.getPeakByMass(mz, new Tolerance(1, true)).setIntensity(1e20f/r);
							
						}
					}
					
					//originalSpec.outputMgf(out2);
					
					spec.outputMgf(out);
				}
				
				if(!train && matlabout != null && spec.getAnnotation() != null){
					
					int i=1;
					for(IonType ion : this.getSigIonsOrderedByIntensityWithOutNoiseIon(spec.getCharge())){
						matlabout.println("MASS(:,"+i+")=[");
						if(ion instanceof IonType.PrefixIon){
							for(float prm : spec.getAnnotation().getPRMMasses(true, 0)){
								matlabout.print(ion.getMz(prm)+"\t");
							}
						}else{
							for(float srm : spec.getAnnotation().getPRMMasses(false, 0)){
								matlabout.print(-ion.getMz(srm)+"\t");
							}
						}
						matlabout.println("];");
						i++;
					}
					
					
					
					matlabout.println("SPEC=[");
					for(Peak p : spec){
						matlabout.println(p.getMz() + "\t" + p.getIntensity());
					}
					matlabout.println("];");
					
					matlabout.println("PROFILE=[");
					for(Peak p : spec){
						if(!profile.containsKey(p)) continue;
						matlabout.print(p.getMz());
						float sum = 0;
						for(IonType ion : this.getSigIonsOrderedByIntensityWithOutNoiseIon(spec.getCharge())){
							float val = Math.max(1e-5f, profile.get(p).get(ion));
							matlabout.print("\t"+ val);
							sum += val;
						}
						matlabout.print("\t"+ (1-sum));
						matlabout.println();
					}
					matlabout.println("];");
				}
				
			}
			out.close();
		}catch (FileNotFoundException e) {
			System.exit(1);
			e.printStackTrace();
		}
		
		System.out.println("Running Done : " + (float)(System.currentTimeMillis() - time)/sn/1000 + " sec/spec");
		
	}

	public void setMaxIterationNum(int maxIterationNum) {
		this.maxIterationNum = maxIterationNum;
	}

	public int getMaxIterationNum() {
		return maxIterationNum;
	}


	public int getMaxCharge() {
		return maxCharge;
	}
	
	public int getMinCharge() {
		return minCharge;
	}


	public static void setForPMCorrection(boolean isForPMCorrection) {
		IPDGenerator.isForPMCorrection = isForPMCorrection;
	}


	public static boolean isForPMCorrection() {
		return isForPMCorrection;
	}
	
}
