package uninovoOld.train;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.IonType;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;
import uninovoOld.IPDGenerator;
import uninovoOld.UniNovo;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.train.InterPeakOffsetFrequencyFunction.InterPeakOffsetPeak;

public class IonSelector {
	private boolean outputIOFF = false;
	final private int minIonNum = 10;
	private float sigprob = 0.15f;
	//private float sigprobForPrecursorIon = 0.2f;
	private int maxIonNum;
	private WindowFilter filter = null;
	private HashMap<IonType, Float> sigIonIntensityMap = null;
	private String specfilename;
	private Tolerance tol;


	public IonSelector(String specfilename, Tolerance tol, WindowFilter filter, int maxIonNum){
		this.specfilename = specfilename;
		this.tol = tol;
		this.filter = filter;
		this.maxIonNum = maxIonNum;
	}
	
	public void setSpecfilename(String specfilename){
		this.specfilename = specfilename;
	}

	private void sortIons(){
		ArrayList<Float> intensities = new ArrayList<Float>();
		HashMap<IonType, Float> tmpMap = new HashMap<IonType, Float>();
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			intensities.add(sigIonIntensityMap.get(ion));
		}
		
		Collections.sort(intensities, Collections.reverseOrder());
		
		for(float intensity : intensities){
		//	if(intensity < 0.15f) continue;
			for(IonType ion : sigIonIntensityMap.keySet()){
				if(intensity == sigIonIntensityMap.get(ion)){
					tmpMap.put(ion, intensity);
					if(tmpMap.size() > maxIonNum-1) break;
				}
			}
			if(tmpMap.size() > maxIonNum-1) break;
		}
		sigIonIntensityMap = tmpMap;
	}
	
	private void normalizeSigIonIntensityMap(){
		
		float maxRatio = 0;
		for(IonType ion : sigIonIntensityMap.keySet()){
			maxRatio = Math.max(maxRatio, sigIonIntensityMap.get(ion));
		}
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			float ratio = sigIonIntensityMap.get(ion);
			ratio /= maxRatio;
			sigIonIntensityMap.put(ion, ratio);
		}
		sortIons();
	}
	
	private void trainIonPeakIntensityRatio(int specCharge){
		Iterator<Spectrum> iterator;
		HashMap<IonType, Float> numMap = new  HashMap<IonType, Float>();
		
		if(sigIonIntensityMap.isEmpty()) return;
		
		int sn = 0;
		ArrayList<IonType> sigIons = new ArrayList<IonType>();
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			if(! (ion instanceof IonType.PrecursorIon)) sigIons.add(ion);
			sigIonIntensityMap.put(ion, 0f);
		}
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 

			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
		
				sn++;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training ion intensities : " + sn);
				}
				
				spec.setRanksOfPeaks();
				PeakGenerator pgen = new PeakGenerator(spec);
				
				
				float highestintensity = 0;
				
				for(Peak p : spec){
					if(p.getRank() == 1){
						highestintensity = p.getIntensity();
						break;
					}
				}
				
				if(filter != null)
					spec = filter.apply(spec);
				
				
				for(Peak p : spec){
					if(p.getIntensity() == 0) continue;
				//	if(p.getRank() > maxRank) continue;
					
					boolean isExplained = false;
					for(IonType ion : sigIons){	
						if(!ion.equals(IonType.NOISE)){
							if(pgen.isExplainedBy(p, ion, tol, tol)){
								Float ratio = sigIonIntensityMap.get(ion);
								Float num = numMap.get(ion);
								
								if(ratio == null) ratio = 0f;
								if(num == null) num = 0f;
								
								ratio += p.getIntensity()/highestintensity;
								num ++;
								
								numMap.put(ion, num);
								sigIonIntensityMap.put(ion, ratio);
								isExplained = true;
							}
						}
					}
					/*if(!isExplained){
						Float ratio = sigIonIntensityMap.get(IonType.NOISE);
						Float num = numMap.get(IonType.NOISE);
						
						if(ratio == null) ratio = 0f;
						if(num == null) num = 0f;
						
						ratio += p.getIntensity()/highestintensity;
						num ++;
						
						numMap.put(IonType.NOISE, num);
						sigIonIntensityMap.put(IonType.NOISE, ratio);
					}*/
				}
			}
			
			for(IonType ion : sigIons){
				float ratio = sigIonIntensityMap.get(ion)/numMap.get(ion);
				sigIonIntensityMap.put(ion, ratio);
			}
			
			normalizeSigIonIntensityMap();
			
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 				
	}
	
	private void trainIonIonCurretRatio(int specCharge){
		Iterator<Spectrum> iterator;
		
		if(sigIonIntensityMap.isEmpty()) return;
		
		int sn = 0;
		ArrayList<IonType> sigIons = new ArrayList<IonType>();
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			if(! (ion instanceof IonType.PrecursorIon)) sigIons.add(ion);
			sigIonIntensityMap.put(ion, 0f);
		}
		
		iterator = UniNovo.getSpectralIterator(specfilename);
		int prevsn = 0; 

		double[] totalIonTic = new double[sigIons.size()];
		double totalTic = 0;
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			
			if(specCharge > 0 && spec.getCharge() != specCharge) continue;
			if(spec.getAnnotation().isModified()) continue;

			sn++;
			if(prevsn < sn/1000){
				prevsn = sn/1000;
				System.out.println("training ion current ratio : " + sn);
			}
		//	if(filter != null) spec  = filter.apply(spec);
			spec.setRanksOfPeaks();
							
			for(int i=0; i<spec.size(); i++){
				Peak p = spec.get(i);
				if(p.getRank() > IPDGenerator.getMaxRank(spec)){
					spec.remove(i--);
				}
			}
			
			
			PeakGenerator pgen = new PeakGenerator(spec);
			
			float tic = 0;
			float[] ionTic = new float[sigIons.size()];
			
			for(Peak p : spec){
				if(p.getIntensity() == 0){
					continue;
				}
				//if(p.getRank() > maxRank){
				//	continue;
				//}
				
				tic += p.getIntensity();
				
				for(IonType ion : sigIons){	
					if(!ion.equals(IonType.NOISE)){
						if(pgen.isExplainedBy(p, ion, tol, tol)){
							ionTic[sigIons.indexOf(ion)] += p.getIntensity();
							continue;
						}
					}
				}
			}
			
			if(tic > 0){
				for(int i=0; i<totalIonTic.length;i++){
					totalIonTic[i] += ionTic[i];
					//System.out.println(totalIonTic[i]);
				}
				totalTic += tic;
			}
		}
		
		for(IonType ion : sigIons){
			sigIonIntensityMap.put(ion, (float)(totalIonTic[sigIons.indexOf(ion)]/totalTic));
		}
		System.out.println(sigIonIntensityMap);
		normalizeSigIonIntensityMap(); 				
	}
	
	private int findSigIons(int specCharge){
		sigIonIntensityMap = new  HashMap<IonType, Float>();
		
		Iterator<Spectrum> iterator;
		HashMap<IonType, Float> tmpSigIonMap = new  HashMap<IonType, Float>();
		HashMap<IonType, Float> sigIonMap = new  HashMap<IonType, Float>();
		
		int sn = 0;
		HashMap<InterPeakOffset, Integer> poffsetsMap = new HashMap<InterPeakOffset, Integer>();
		HashMap<InterPeakOffset, Integer> soffsetsMap = new HashMap<InterPeakOffset, Integer>();
		HashMap<InterPeakOffset, Integer> roffsetsMap = new HashMap<InterPeakOffset, Integer>();
		
		iterator = UniNovo.getSpectralIterator(specfilename);
		int prevsn = 0; 
		float normalizer = 0;
		float normalizerForPrecursor = 0;
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			
			spec.setRanksOfPeaks();
			
			for(int i=0; i<spec.size(); i++){
				Peak p = spec.get(i);
				if(p.getRank() > IPDGenerator.getMaxRank(spec)){
					spec.remove(i--);
				}
			}
			
			
			Spectrum filteredSpec = spec;
			if(filter != null)filteredSpec = filter.apply(spec);
			if(specCharge > 0 && spec.getCharge() != specCharge) continue;
			if(spec.getAnnotation().isModified()) continue;

			sn++;
			if(prevsn < sn/1000){
				prevsn = sn/1000;
				System.out.println("finding significant ions : " + sn);
			}
			
			spec.setRanksOfPeaks();
			PeakGenerator pgen = new PeakGenerator(spec);
			
			
			for(int charge = 1; charge <= spec.getCharge(); charge++){
				
				ArrayList<Peak> comparedPeaks = new ArrayList<Peak>();
				float maxMz = PeakParameter.maxMzWith(charge, spec);
				for(Peak p : filteredSpec){
					if(p.getIntensity() > 0 && p.getMz() < maxMz){
						comparedPeaks.add(p);
					}
				}
				
				
				ArrayList<Peak> pbps = pgen.getTheoreticalPrefixBasePeaks(charge);
				ArrayList<Peak> sbps = pgen.getTheoreticalSuffixBasePeaks(charge);
				Peak rbp = pgen.getTheoreticalPrecursorBasePeak(charge);
				
				if(charge == 1){
					normalizer += pbps.size();
					normalizerForPrecursor ++;
				}
				
				ArrayList<InterPeakOffset> poffsets = new ArrayList<InterPeakOffset>();
				ArrayList<InterPeakOffset> soffsets = new ArrayList<InterPeakOffset>();
				ArrayList<InterPeakOffset> roffsets = new ArrayList<InterPeakOffset>();
				
				for(Peak bp : pbps)
					poffsets.addAll(InterPeakOffsetFrequencyFunction.getGeneralizedOffsets(bp, comparedPeaks, false, 0, tol));
				for(Peak bp : sbps)
					soffsets.addAll(InterPeakOffsetFrequencyFunction.getGeneralizedOffsets(bp, comparedPeaks, false, 0, tol));
				roffsets.addAll(InterPeakOffsetFrequencyFunction.getGeneralizedOffsets(rbp, comparedPeaks, false, 0, tol));
					
			/*	if(poffsets.size() > pbps.size()){
					System.out.println("**" + poffsets.size() + "\t" + pbps.size() + "\t" + spec.getAnnotationStr());
					for(InterPeakOffset poff : poffsets)
						System.out.println(poff);
				}*/
				updateOffsetsMap(poffsetsMap, poffsets);
				updateOffsetsMap(soffsetsMap, soffsets);
				updateOffsetsMap(roffsetsMap, roffsets);//TODO something is wrong with precursor ion selection
				
				
			}
		}

		String iofffilename = null;
		if(outputIOFF) iofffilename = specfilename.substring(0, specfilename.length()-4) + "_"  + specCharge + "_" + "Prefix";
		
		for(InterPeakOffsetPeak gofPeak : InterPeakOffsetFrequencyFunction.getOffSetFrequencyFunction(poffsetsMap, normalizer, 0 , iofffilename)){
			InterPeakOffset gof = gofPeak.getInterPeakOffset();
			tmpSigIonMap.put(new IonType.PrefixIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getProbability());
		}	

		if(outputIOFF) iofffilename = specfilename.substring(0, specfilename.length()-4) + "_"  + specCharge + "_"  + "Suffix";
		for(InterPeakOffsetPeak gofPeak : InterPeakOffsetFrequencyFunction.getOffSetFrequencyFunction(soffsetsMap, normalizer, 0, iofffilename)){
			InterPeakOffset gof = gofPeak.getInterPeakOffset();
			tmpSigIonMap.put(new IonType.SuffixIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getProbability());
		}	
		
		if(outputIOFF) iofffilename = specfilename.substring(0, specfilename.length()-4) + "_"  + specCharge + "_"  + "Precursor";
		for(InterPeakOffsetPeak gofPeak : InterPeakOffsetFrequencyFunction.getOffSetFrequencyFunction(roffsetsMap, normalizerForPrecursor, 0, iofffilename)){
			InterPeakOffset gof = gofPeak.getInterPeakOffset();
			tmpSigIonMap.put(new IonType.PrecursorIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getProbability());
		}	
		
		
		ArrayList<IonType> knownIonTypes = IonType.getAllKnownIonTypes(Math.min(4,specCharge), false);
		
		for(IonType ion : tmpSigIonMap.keySet()){
			if(ion instanceof IonType.PrecursorIon) continue;
			
			boolean isKnown = false;
			float diff = 10000f;
			
			IonType keyIon = null;
			float v = 0;
			for(IonType kion : knownIonTypes){
				if(ion.getCharge() != kion.getCharge()) continue;
				if(kion instanceof IonType.PrecursorIon) continue;
				
				if((ion.isPrefixIon() && kion.isPrefixIon()) || (!ion.isPrefixIon() && !kion.isPrefixIon())){
					float tdiff = Math.abs(ion.getOffset() - kion.getOffset());
					
					if(tdiff < diff){
						diff = tdiff;
						if(diff < 0.5){
							Float prev = sigIonMap.get(kion);
							if(prev == null) prev = 0f;
							
							keyIon = kion;
							v = Math.max(prev, tmpSigIonMap.get(ion));
						
							isKnown = true;
						}
					}
				}
			}
			if(!isKnown) sigIonMap.put(ion, tmpSigIonMap.get(ion));
			else sigIonMap.put(keyIon, v);
			
		}

		if(!sigIonMap.isEmpty()){
			ArrayList<Float> probs = new ArrayList<Float>();
			for(IonType ion : sigIonMap.keySet()){
				probs.add(sigIonMap.get(ion));
			}
			
			Collections.sort(probs);
			
			for(IonType ion : sigIonMap.keySet()){
				if(ion instanceof IonType.PrecursorIon){
				//	if(sigIonMap.get(ion) >= sigprobForPrecursorIon)
				//		sigIonIntensityMap.put(ion, 0f);
				}else	if(sigIonMap.get(ion) >= Math.min(sigprob, probs.get(Math.max(0, probs.size()-minIonNum))))
				//	System.out.println(ion + "\t" + sigIonMap.get(ion));
					sigIonIntensityMap.put(ion,sigIonMap.get(ion));
			}
		//	sigIonIntensityMap.put(IonType.NOISE, 0f);
			//System.out.println(sigIonIntensityMap);
			
			sortIons();
			
			print(sigIonIntensityMap);
			
			normalizeSigIonIntensityMap();
			//System.out.println(sigIonIntensityMap);
		} 
		return sn;
	}
	
	static private void print(HashMap<IonType, Float> ionMap){
		ArrayList<Float> f = new ArrayList<Float>(ionMap.values());
		Collections.sort(f, Collections.reverseOrder());
		
		System.out.print("\\multirow{2}{*}{CID2}&$\\delta$&");
		for(float i : f){
			for(IonType ion : ionMap.keySet()){
				if(ionMap.get(ion) == i){
					System.out.print("$"+ion.getName() + "$&");
				}
			}
		}
		System.out.println();
		System.out.print("&$OFF(\\delta)$&");
		for(float i : f){
			for(IonType ion : ionMap.keySet()){
				if(ionMap.get(ion) == i){
					System.out.print(String.format("%.2f", i)  + "&");
				}
			}
		}
		System.out.println();
	}
	
	static private void updateOffsetsMap(HashMap<InterPeakOffset, Integer> offsetsMap, ArrayList<InterPeakOffset> offsets){
		if(offsets==null) return;

		for(InterPeakOffset gof : offsets){
			Integer n = offsetsMap.get(gof);
			if(n==null) n = 0;
			
			n++;
			
			offsetsMap.put(gof, n);
		}
	}
	
	public int train(int specCharge, boolean considerIonIntensity){
		int sn = findSigIons(specCharge);
		if(considerIonIntensity) trainIonIonCurretRatio(specCharge);//trainIonPeakIntensityRatio(specCharge);
		
		/*int n = 0;
		for(IonType ion : sigIonIntensityMap.keySet()){
			if(ion instanceof IonType.PrecursorIon) continue;
			n++;

		}
		*/
		//print(sigIonIntensityMap);
		//System.out.println(sigIonIntensityMap);
		//System.out.println(n);
		
		return sn;
	}
	
	public HashMap<IonType, Float> getSigIons(){
		return sigIonIntensityMap;
	}

	static public void main(String[] args){
		String inputmgf = "/home/kwj/workspace/inputs/Training/ETDAspN_train.mgf";
		int charge = 4;
		int maxIonNum = 8;
		Tolerance tol = new Tolerance(0.5f, false);
		//tol = new Tolerance(20f, true);
		IonSelector  ionfinder = new IonSelector(inputmgf, tol, new WindowFilter(6, 50), maxIonNum);
		ionfinder.train(charge, false);
		
		
	}
	
}
