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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import uninovo.FeatureFrequencyFunction.FeatureFrequencyFunctionPeak;
import uninovo.util.AminoAcidSet;
import uninovo.util.IonType;
import uninovo.util.Peak;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
import uninovo.util.WindowFilter;

/**
 * The Class FeatureSelector selects significant features in the training session.
 */
public class FeatureSelector {	
	
	/** The feature map. */
	private HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> featureMap = null;
	
	/** The window filter to filter out peaks in spectra. */
	private WindowFilter filter = null;
	
	/** The MS1 tol. */
	private Tolerance pmtol;
	
	/** The significant ions. */
	private ArrayList<IonType> sigIons;
	
	/** The probability threshold to select a feature. */
	public float sigprob = 0.10f;
	
	/** The specfile name. */
	private String specfilename;
	
	/** The MS2 tol. */
	private Tolerance tol;
	
	/**
	 * Instantiates a new feature selector.
	 *
	 * @param specfilename the spec filename
	 * @param sigIons the significant ions
	 * @param tol MS2 tol
	 * @param pmtol MS1 tol
	 */
	public FeatureSelector(String specfilename, ArrayList<IonType> sigIons, Tolerance tol, Tolerance pmtol){
		this.specfilename = specfilename;
		this.sigIons = sigIons;
		this.tol = tol;
		this.pmtol = pmtol;
	}
	
	/**
	 * Filter.
	 *
	 * @param b window filter b
	 * @return this
	 */
	public FeatureSelector filter(WindowFilter b) {filter = b; return this;}
	
	/**
	 * Find significant features
	 *
	 * @param sigFeatures the significant features
	 * @param ion the ion
	 * @param specCharge the spec charge
	 * @param iterationNum the iteration num
	 */
	private void  
		findSigFeatures(HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<OffsetFeature>>> 
		sigFeatures,
		IonType ion,
		int specCharge, int iterationNum){
		
		Iterator<Spectrum> iterator;
		HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>>> 
			ffMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>>>();
		HashMap<SpectrumParameter, HashMap<PeakParameter, Integer>> 
			numMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, Integer>>();
	
		iterator = UniNovo.getSpectralIterator(specfilename);
		int prevsn = 0; int sn = 0;
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			Spectrum filteredspec = spec;
			int maxRank = PeakIonProbabilityGenerator.getMaxRank(spec);
			
			if(spec.getCharge() != specCharge) continue;
			if(spec.getAnnotation().isModified()) continue;
			
			
			sn++;
			if(prevsn < sn/1000){
				prevsn = sn/1000;
				System.out.println("Iteration: " +iterationNum+ " selecting significant features for " + ion + ": " + sn);
			}
			
			spec.setRanksOfPeaks();
			
			for(int i=0; i<spec.size(); i++){
				Peak p = spec.get(i);
				if(p.getRank() > maxRank){
					spec.remove(i--);
				}
			}
			
			if(filter != null){
				filteredspec = filter.apply(spec);
			}
			
			
			SpectrumParameter spar = new SpectrumParameter(spec);
			PeakGenerator pgen = new PeakGenerator(spec);
			
			if(!numMap.containsKey(spar)) numMap.put(spar,  new HashMap<PeakParameter, Integer>());
			HashMap<PeakParameter, Integer> nums = numMap.get(spar);
			
			for(Peak bp : spec){
				if(bp.getRank() > maxRank) continue;
					
				if(ion.equals(IonType.NOISE)){
					boolean expd = false;
					for(IonType i : sigIons) if(pgen.isExplainedBy(bp, i, tol, pmtol)) {expd = true; break;}
					if(expd) continue;
				}else	if(!pgen.isExplainedBy(bp, ion, tol, pmtol)) continue;
				
				PeakParameter ppar = new PeakParameter(bp, spec, iterationNum);
				Integer n = nums.get(ppar);
				if(n == null) n = 0;
				nums.put(ppar, n+1);
				
				int charge = ion.getCharge();
				
				for(int chargeOffset=1 - charge; chargeOffset<=spec.getCharge() - charge; chargeOffset++){
					for(int i=0;i<2;i++){
						Peak nbp = PeakGenerator.getChargeChangedPeak(bp, charge, chargeOffset);
						if(i==1) nbp = PeakGenerator.getComplementaryPeak(nbp, nbp.getCharge(), spec);
						
						float minMass = nbp.getMz() + FeatureFrequencyFunction.MIN/nbp.getCharge();
						float maxMass = nbp.getMz() + FeatureFrequencyFunction.MAX/nbp.getCharge();
						
						ArrayList<Peak> compPeaks = filteredspec.getPeakListByMassRange(minMass, maxMass);
						if(nbp.getCharge() == filteredspec.getCharge())
							compPeaks.add(filteredspec.getPrecursorPeak());
						
						for(Peak cp : compPeaks){ 
							if(cp.equals(bp)) continue;			 
						//	if(cp.getRank() > maxRank) continue;
							int peakIntensityRatio = PeakParameter.getPeakIntensityRatioNum(bp, cp, spec);

							RelationBetweenPeaks off = FeatureFrequencyFunction.getRelationsBetween(nbp, cp, i==1, chargeOffset, tol);
							updateFFFMap(ffMap, spar, ppar, peakIntensityRatio, off);
						}
					}
				}
			}	
		}
		
		updateSigFeatures(sigFeatures, ffMap, numMap, iterationNum, ion); 
	}
	
	/**
	 * Gets the significant features.
	 *
	 * @return the significant features
	 */
	public HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> getSignificantFeatures(){
		return featureMap;
	}

	
	/**
	 * Train to select significant features
	 *
	 * @param specCharge the spectrum charge
	 * @param iterationNum the iteration number
	 * @param aaSet the amino acid set
	 */
	public void train(int specCharge, int iterationNum, AminoAcidSet aaSet){
		featureMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>>();
		HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<OffsetFeature>>>
		sigFeatures = new HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<OffsetFeature>>>();
			
		for(IonType ion : sigIons){
			if(ion instanceof IonType.PrecursorIon) continue;
			findSigFeatures(sigFeatures, ion, specCharge, iterationNum);
		}
		
		//findSigGOFCons(sigGofCons, IonType.NOISE, specCharge, 30, iterationNum);
		
		HashSet<Integer> charges = new HashSet<Integer>();
		
		for(IonType ion : sigIons){
			charges.add(ion.getCharge());
		}
		
		for(SpectrumParameter spar : sigFeatures.keySet()){
			if(!featureMap.containsKey(spar)) featureMap.put(spar, new HashMap<PeakParameter, ArrayList<Feature>>());
			HashMap<PeakParameter, ArrayList<Feature>> subFeatureMap = featureMap.get(spar);
		
			for(PeakParameter ppar : sigFeatures.get(spar).keySet()){
				if(!subFeatureMap.containsKey(ppar)) subFeatureMap.put(ppar, new ArrayList<Feature>());
				ArrayList<Feature> features = subFeatureMap.get(ppar);
				
				features.addAll(sigFeatures.get(spar).get(ppar));
			}
		}
		
		float[] n = new float[2];
		for(SpectrumParameter spar : SpectrumParameter.getAllSpectrumParameters(specCharge)){
			if(!featureMap.containsKey(spar)) featureMap.put(spar, new HashMap<PeakParameter, ArrayList<Feature>>());
			HashMap<PeakParameter, ArrayList<Feature>> subFeatureMap = featureMap.get(spar);
			
			for(PeakParameter ppar : PeakParameter.getAllBasePeakParameters(spar.getSpecCharge(), iterationNum)){
				if(!subFeatureMap.containsKey(ppar)) subFeatureMap.put(ppar, new ArrayList<Feature>());
				ArrayList<Feature> features = subFeatureMap.get(ppar);
				for(int charge : charges){
					for(int ratio : PeakParameter.getAllPeakIntensityRatioNums()){
						features.add(new LinkingFeature(spar, ppar, charge, true, iterationNum, ratio, aaSet));
						features.add(new LinkingFeature(spar, ppar, charge, false, iterationNum, ratio, aaSet));
					
					}
					features.add(new IntensityFeature(spar, ppar, charge, iterationNum));
				}
				n[0] ++; n[1] += features.size();
			}
		}
		System.out.println("Average # features per one peak: " + n[1] / n[0]);
	}
	
	/**
	 * Update ff (feature offset) map.
	 *
	 * @param ffMap the ff map
	 * @param spar the spec param
	 * @param ppar the peak param
	 * @param peakIntensityRatio the peak intensity ratio
	 * @param ff the ff
	 */
	private void updateFFFMap(HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>>> 
	 ffMap, SpectrumParameter spar, PeakParameter ppar, int peakIntensityRatio, RelationBetweenPeaks ff){
		
		if(ff == null) return;
		
		if(!ffMap.containsKey(spar)) ffMap.put(spar, new HashMap<PeakParameter, HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>>());
		HashMap<PeakParameter, HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>> subFFMap1 = ffMap.get(spar);
	
		if(!subFFMap1.containsKey(ppar)) subFFMap1.put(ppar,  new HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>());
		HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>> subFFMap2 = subFFMap1.get(ppar);
	
		if(!subFFMap2.containsKey(peakIntensityRatio)) subFFMap2.put(peakIntensityRatio, new HashMap<RelationBetweenPeaks, Integer>());
		HashMap<RelationBetweenPeaks, Integer> subFFMap3 = subFFMap2.get(peakIntensityRatio);
		
		Integer num = subFFMap3.get(ff);
		if(num == null) num = 0;
		num ++;
		
		subFFMap3.put(ff, num);
	}
	
	/**
	 * Update significant features
	 *
	 * @param sigFeatures the sig features
	 * @param ffMap the ff map
	 * @param numMap the num map
	 * @param iterationNum the iteration num
	 * @param ion the ion
	 */
	private void updateSigFeatures(
		HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<OffsetFeature>>>
		sigFeatures,
		HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>>> 
		ffMap, 
		HashMap<SpectrumParameter, HashMap<PeakParameter, Integer>> 
		numMap,
		int iterationNum,
		IonType ion)
	{
		for(SpectrumParameter spar : ffMap.keySet()){
			HashMap<PeakParameter, HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>>> subGofMap1 = ffMap.get(spar);
			if(!sigFeatures.containsKey(spar)) sigFeatures.put(spar, new HashMap<PeakParameter, ArrayList<OffsetFeature>>());
			HashMap<PeakParameter, ArrayList<OffsetFeature>> subSigGofCons = sigFeatures.get(spar);
		
			for(PeakParameter ppar : subGofMap1.keySet()){
				HashMap<Integer, HashMap<RelationBetweenPeaks, Integer>> subGofMap2 = subGofMap1.get(ppar);
				if(!subSigGofCons.containsKey(ppar)) subSigGofCons.put(ppar, new ArrayList<OffsetFeature>());
				ArrayList<OffsetFeature> gofcs = subSigGofCons.get(ppar);
			
				for(int peakIntensityRatio : subGofMap2.keySet()){			
					for(FeatureFrequencyFunctionPeak gofPeak : FeatureFrequencyFunction.getFeatureFrequencyFunction(subGofMap2.get(peakIntensityRatio), numMap.get(spar).get(ppar), sigprob)){
						RelationBetweenPeaks gof = gofPeak.getFFF();
						OffsetFeature gofConp = new OffsetFeature(spar, ppar, gof.getBaseCharge(), true, iterationNum, peakIntensityRatio, gof);
						OffsetFeature gofCona = new OffsetFeature(spar, ppar, gof.getBaseCharge(), false, iterationNum, peakIntensityRatio, gof);
						
						if(!gofcs.contains(gofConp)) gofcs.add(gofConp);
						if(!gofcs.contains(gofCona)) gofcs.add(gofCona);
					}
				}
			}
		}
	}
}
