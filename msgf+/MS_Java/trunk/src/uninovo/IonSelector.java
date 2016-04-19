/***************************************************************************
  * Title:          
  * Author:         Kyowon Jeong
  * Last modified:  
  *
  * Copyright (c) 2008-2009 The Regents of the University of California
  * All Rights Reserved
  * See file LICENSE for details.
  ***************************************************************************/
package uninovo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import uninovo.FeatureFrequencyFunction.FeatureFrequencyFunctionPeak;
import uninovo.util.IonType;
import uninovo.util.Peak;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
import uninovo.util.WindowFilter;

/**
 * The Class IonSelector selects significant ion types using offset frequency function (OFF)
 */
public class IonSelector {
	
	/**
	 * Update feature map.
	 *
	 * @param featureMap the feature map
	 * @param features the features
	 */
	static private void updateFeatureMap(HashMap<RelationBetweenPeaks, Integer> featureMap, ArrayList<RelationBetweenPeaks> features){
		if(features==null) return;

		for(RelationBetweenPeaks gof : features){
			Integer n = featureMap.get(gof);
			if(n==null) n = 0;
			
			n++;
			
			featureMap.put(gof, n);
		}
	}
	
	/** The filter. */
	private WindowFilter filter = null;
	
	/** The max ion num. */
	private int maxIonNum;
	
	/** The min ion num. */
	final private int minIonNum = 10;
	
	/** The significant ion intensity map. */
	private HashMap<IonType, Float> sigIonIntensityMap = null;
	
	/** The significant prob to select ions. */
	private float sigprob = 0.15f;
	
	/** The spec filename. */
	private String specfilename;


	/** The MS2 tol. */
	private Tolerance tol;

	/**
	 * Instantiates a new ion selector.
	 *
	 * @param specfilename the spec filename
	 * @param tol the MS1 tol
	 * @param filter the window filter
	 * @param maxIonNum the max ion number
	 */
	public IonSelector(String specfilename, Tolerance tol, WindowFilter filter, int maxIonNum){
		this.specfilename = specfilename;
		this.tol = tol;
		this.filter = filter;
		this.maxIonNum = maxIonNum;
	}
	
	/**
	 * Find significant ions.
	 *
	 * @param specCharge the spec charge
	 * @return the number of spectra
	 */
	private int findSigIons(int specCharge){
		sigIonIntensityMap = new  HashMap<IonType, Float>();
		
		Iterator<Spectrum> iterator;
		HashMap<IonType, Float> tmpSigIonMap = new  HashMap<IonType, Float>();
		HashMap<IonType, Float> sigIonMap = new  HashMap<IonType, Float>();
		
		int sn = 0;
		HashMap<RelationBetweenPeaks, Integer> poffsetsMap = new HashMap<RelationBetweenPeaks, Integer>();
		HashMap<RelationBetweenPeaks, Integer> soffsetsMap = new HashMap<RelationBetweenPeaks, Integer>();
		HashMap<RelationBetweenPeaks, Integer> roffsetsMap = new HashMap<RelationBetweenPeaks, Integer>();
		
		iterator = UniNovo.getSpectralIterator(specfilename);
		int prevsn = 0; 
		float normalizer = 0;
		float normalizerForPrecursor = 0;
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			
			spec.setRanksOfPeaks();
			
			for(int i=0; i<spec.size(); i++){
				Peak p = spec.get(i);
				if(p.getRank() > PeakIonProbabilityGenerator.getMaxRank(spec)){
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
				
				
				ArrayList<Peak> pbps = pgen.getTheoreticalPrefixPeaks(charge);
				ArrayList<Peak> sbps = pgen.getTheoreticalSuffixPeaks(charge);
				Peak rbp = pgen.getTheoreticalPrecursorPeak(charge);
				
				if(charge == 1){
					normalizer += pbps.size();
					normalizerForPrecursor ++;
				}
				
				ArrayList<RelationBetweenPeaks> poffsets = new ArrayList<RelationBetweenPeaks>();
				ArrayList<RelationBetweenPeaks> soffsets = new ArrayList<RelationBetweenPeaks>();
				ArrayList<RelationBetweenPeaks> roffsets = new ArrayList<RelationBetweenPeaks>();
				
				for(Peak bp : pbps)
					poffsets.addAll(FeatureFrequencyFunction.getFeaturesBetween(bp, comparedPeaks, false, 0, tol));
				for(Peak bp : sbps)
					soffsets.addAll(FeatureFrequencyFunction.getFeaturesBetween(bp, comparedPeaks, false, 0, tol));
				roffsets.addAll(FeatureFrequencyFunction.getFeaturesBetween(rbp, comparedPeaks, false, 0, tol));
					
	
				updateFeatureMap(poffsetsMap, poffsets);
				updateFeatureMap(soffsetsMap, soffsets);
				updateFeatureMap(roffsetsMap, roffsets);
				
			}
		}
		
		for(FeatureFrequencyFunctionPeak gofPeak : FeatureFrequencyFunction.getFeatureFrequencyFunction(poffsetsMap, normalizer, 0)){
			RelationBetweenPeaks gof = gofPeak.getFFF();
			tmpSigIonMap.put(new IonType.PrefixIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getFrequency());
		}	

		for(FeatureFrequencyFunctionPeak gofPeak : FeatureFrequencyFunction.getFeatureFrequencyFunction(soffsetsMap, normalizer, 0)){
			RelationBetweenPeaks gof = gofPeak.getFFF();
			tmpSigIonMap.put(new IonType.SuffixIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getFrequency());
		}	
		
		for(FeatureFrequencyFunctionPeak gofPeak : FeatureFrequencyFunction.getFeatureFrequencyFunction(roffsetsMap, normalizerForPrecursor, 0)){
			RelationBetweenPeaks gof = gofPeak.getFFF();
			tmpSigIonMap.put(new IonType.PrecursorIon(gof.getBaseCharge(), gof.getOffset()/gof.getBaseCharge()), gofPeak.getFrequency());
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
				}else	if(sigIonMap.get(ion) >= Math.min(sigprob, probs.get(Math.max(0, probs.size()-minIonNum))))
					sigIonIntensityMap.put(ion,sigIonMap.get(ion));
			}
			
			sortIons();
			
			//print(sigIonIntensityMap);
			
			normalizeSigIonIntensityMap();
		} 
		return sn;
	}
	
	/**
	 * Gets the sig ions.
	 *
	 * @return the significant ions
	 */
	public HashMap<IonType, Float> getSigIons(){
		return sigIonIntensityMap;
	}
	
	/**
	 * Normalize significant ion intensity map.
	 */
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
	
	/**
	 * Sort ions.
	 */
	private void sortIons(){
		ArrayList<Float> intensities = new ArrayList<Float>();
		HashMap<IonType, Float> tmpMap = new HashMap<IonType, Float>();
		
		for(IonType ion : sigIonIntensityMap.keySet()){
			intensities.add(sigIonIntensityMap.get(ion));
		}
		
		Collections.sort(intensities, Collections.reverseOrder());
		
		for(float intensity : intensities){
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
	
	/**
	 * Train.
	 *
	 * @param specCharge the spec charge
	 * @param considerIonIntensity the consider ion intensity
	 * @return the number of spectra
	 */
	public int train(int specCharge, boolean considerIonIntensity){
		int sn = findSigIons(specCharge);
		if(considerIonIntensity) trainIonCurrentRatio(specCharge);
		return sn;
	}
	
	/**
	 * Train ion current ratio.
	 *
	 * @param specCharge the spec charge
	 */
	private void trainIonCurrentRatio(int specCharge){
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
			spec.setRanksOfPeaks();
							
			for(int i=0; i<spec.size(); i++){
				Peak p = spec.get(i);
				if(p.getRank() > PeakIonProbabilityGenerator.getMaxRank(spec)){
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

	
}
