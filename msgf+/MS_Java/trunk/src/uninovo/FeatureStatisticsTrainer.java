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

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import uninovo.util.IonType;
import uninovo.util.Peak;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
import uninovo.util.WindowFilter;

/**
 * The Class FeatureStatisticsTrainer learns probabilities associated with features
 */
public class FeatureStatisticsTrainer {
	
	/** The KL dthreshold. */
	static float KLDthreshold = 0;
	
	/** The min feature num per group. */
	static int minFeatureNumPerGroup = 50;
	
	/** The feature map. */
	private HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> featureMap;
	
	/** The filter. */
	private WindowFilter filter = null;
	
	/** The independent feature map. */
	private HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>> independentFeatureMap;
	
	/** The MS1 tol. */
	private Tolerance pmtol;
	
	/** The significant ion map. */
	private HashMap<IonType, Float> sigIonMap;
	
	/** The significant ions. */
	private ArrayList<IonType> sigIons;
	
	/** The spec filename. */
	private String specfilename;
	
	/** The MS2 tol. */
	private Tolerance tol;
	
	/**
	 * Instantiates a new feature statistics trainer.
	 *
	 * @param specfilename the spec filename
	 * @param sigIonMap the sig ion map
	 * @param featureMap the feature map
	 * @param tol the MS2 tol
	 * @param pmtol the MS1 tol
	 */
	public FeatureStatisticsTrainer(String specfilename,  HashMap<IonType, Float> sigIonMap, HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> featureMap,
			Tolerance tol, Tolerance pmtol){
		this.specfilename = specfilename;
		this.sigIonMap = sigIonMap;
		sigIons = new ArrayList<IonType>();
		
		ArrayList<Float> intensities = new ArrayList<Float>();
		
		for(IonType ion : sigIonMap.keySet()){
			intensities.add(sigIonMap.get(ion));
		}
		
		Collections.sort(intensities);
		
		for(int i = intensities.size()-1; i>=0; i--){
			for(IonType ion : sigIonMap.keySet()){
				if(intensities.get(i) == sigIonMap.get(ion) && !ion.equals(IonType.NOISE)){
					sigIons.add(ion);
				}
			}
		}
		
		this.featureMap = featureMap;
		this.tol = tol;
		this.pmtol = pmtol;
	}
	
	/**
	 * Discard features with low divergences
	 *
	 * @param minCaseNum the min case num
	 * @param maxFeatureNum the max feature num
	 * @param discard to discard or not
	 */
	public void discardFeaturesWithLowDivergences(int minCaseNum, int maxFeatureNum, boolean discard){
		ArrayList<Feature> toErase = new ArrayList<Feature>();
		
		int[] num = new int[4];
		int[] erasednum = new int [3];
		ArrayList<Float> KLDs = new ArrayList<Float>();
		
		for(SpectrumParameter spar : featureMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> cm = featureMap.get(spar);
			for(PeakParameter ppar : cm.keySet()){
				for(Feature feature : cm.get(ppar)){
					if(feature instanceof IntensityFeature){
						num[0]++;
						feature.calculateProbabilities();
					}
				}
			}
		}
	
		for(SpectrumParameter spar : featureMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> subFeatureMap = featureMap.get(spar);
			for(PeakParameter ppar : subFeatureMap.keySet()){
				ArrayList<Feature> features = subFeatureMap.get(ppar);
				for(int i=0; i<features.size();i++){
					Feature feature = features.get(i);
					if(feature instanceof IntensityFeature) continue;
					else if(feature instanceof OffsetFeature) num[1]++;
					else if(feature instanceof LinkingFeature) num[2]++;
					
					float sum = feature.calculateProbabilities();
					if(discard && sum < minCaseNum){
						features.remove(i--);
						toErase.add(feature);
					}else{
						KLDs.add(feature.getDivergence());
					}
				}
			}
		}
		
		Collections.sort(KLDs, Collections.reverseOrder());
		KLDthreshold = 0;

		if(!KLDs.isEmpty())
			KLDthreshold = KLDs.get(Math.min(KLDs.size()-1, maxFeatureNum));
		
		for(SpectrumParameter spar : featureMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> subFeatureMap = featureMap.get(spar);
			for(PeakParameter ppar : subFeatureMap.keySet()){
				ArrayList<Feature> features = subFeatureMap.get(ppar);
				for(int i=0; i<features.size();i++){
					Feature feature = features.get(i);
					if(feature instanceof IntensityFeature) continue;
				
					if(discard && feature.getDivergence() < KLDthreshold){					
					}
			
				}
			}
		}
		
		
		for(Feature con : toErase){
			if(con instanceof OffsetFeature) erasednum[0]++;
			else if(con instanceof LinkingFeature) erasednum[1]++;	
		}
		
		System.out.println("# Null feature : " + num[0]);
		System.out.println("# Ion dependency feature : " + num[1] + " -> " + (num[1] - erasednum[0]));
		System.out.println("# Linking feature : " + num[2] + " -> " + (num[2] - erasednum[1]));

	}
	
	/**
	 * Filter.
	 *
	 * @param b the window filter
	 * @return this
	 */
	public FeatureStatisticsTrainer filter(WindowFilter b) {filter = b; return this;}
	
	
	/**
	 * Train the statistics
	 *
	 * @param specCharge the spec charge
	 * @param iterationNum the iteration number
	 */
	public void train(int specCharge,  int iterationNum){
		
		int sn = 0;
		Iterator<Spectrum> iterator = UniNovo.getSpectralIterator(specfilename);
		int prevsn = 0;
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			Spectrum filteredspec = spec;
			
			if(spec.getCharge() != specCharge) continue;
			if(spec.getAnnotation().isModified()) continue;
			
			
			sn++;
			
			if(prevsn < sn/1000){
				prevsn = sn/1000;
				System.out.println("Iteration: " +iterationNum+ " Feature probability training " +": " + sn);
			}
			
			spec.setRanksOfPeaks();
			int maxRank = PeakIonProbabilityGenerator.getMaxRank(spec);
			
			for(int i=0; i<spec.size(); i++){
				Peak p = spec.get(i);
				if(p.getRank() > maxRank){
					spec.remove(i--);
				}
			}
			

			if(filter!=null){
				filteredspec = filter.apply(spec);
			}
			
			PeakGenerator pgen = new PeakGenerator(spec);
			SpectrumParameter spar = new SpectrumParameter(spec);	
			HashMap<PeakParameter, ArrayList<Feature>> sigFeatureMap = featureMap.get(spar);
			
			for(Peak bp : spec){
				if(bp.getRank() > maxRank) continue;
				
				PeakParameter ppar = new PeakParameter(bp, spec,iterationNum);
				ArrayList<Feature> sigFeatures = sigFeatureMap.get(ppar);
				ArrayList<Feature> matchedFeatures = new ArrayList<Feature>();
				IntensityFeature intensityfeature = null;
				IonType expIon = null;
				
				for(IonType ion : sigIons){
					if(pgen.isExplainedBy(bp, ion, tol, pmtol)){
						expIon = ion;
						break;
					}
				}
				
				if(expIon == null) expIon = IonType.NOISE;
				
				for(Feature feature : sigFeatures){
					if(feature.isSatisfiedBy(bp, filteredspec, spar, ppar, tol, pmtol, iterationNum)){
						matchedFeatures.add(feature);
					}
					if(feature instanceof IntensityFeature){
						intensityfeature = (IntensityFeature) feature;			
					}	
				}
				
				for(Feature feature : matchedFeatures){
					feature.setIntensityFeature(intensityfeature);
					feature.addIonCount(expIon);
				}
			}
		}
	}
	
	/**
	 * Update independent feature groups.
	 */
	private void updateIndependentFeatureGroups(){
		independentFeatureMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>();

		for(SpectrumParameter spar : featureMap.keySet()){
			int sc = spar.getSpecCharge();
			
			HashMap<PeakParameter, ArrayList<Feature>> subFeatureMap = featureMap.get(spar);
			
			if(!independentFeatureMap.containsKey(spar)) independentFeatureMap.put(spar,  new HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>());
			HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> si = independentFeatureMap.get(spar);
				
			for(PeakParameter ppar : subFeatureMap.keySet()){
				if(!si.containsKey(ppar)) si.put(ppar, new HashMap<Integer, ArrayList<Feature>>());
				HashMap<Integer, ArrayList<Feature>> ssi = si.get(ppar);
				for(Feature feature : subFeatureMap.get(ppar)){
					int id;
					if(feature instanceof IntensityFeature) id = 0;
					else if(feature instanceof LinkingFeature){
						id = feature.getBasePeakCharge();
					}else{
						OffsetFeature gcon = (OffsetFeature) feature;
											
						if(!gcon.isComplementary() && gcon.getChargeOffset() == 0) id = sc + 1;
						else{
							id = sc + 2 + (gcon.isComplementary()? sc * sc : 0) + sc * (gcon.getBasePeakCharge() - 1) + (gcon.getBasePeakCharge() + gcon.getChargeOffset() - 1);
						}
					}
					if(!ssi.containsKey(id)) ssi.put(id, new ArrayList<Feature>());
					ArrayList<Feature> features = ssi.get(id);
					
					features.add(feature);
				}
			}	
		}
	}
	
	
	/**
	 * Write the parameters in file.
	 *
	 * @param file the file
	 * @param charge the charge
	 * @param iterationNum the iteration number
	 */
	public void writeInFile(String file, int charge, int iterationNum){
		updateIndependentFeatureGroups();
		
		PrintWriter out;
		try {
			out = new PrintWriter(new FileWriter(file, true));
			
			out.println("#SPECPARTITION\t" + charge);
			out.println(SpectrumParameter.writePartitionMzs(charge));
			
			out.println("#PMTOL\t" + pmtol.toString());
			out.println("#TOL\t" + tol.toString());
			out.println("#MAXSPECMZ\t"+charge +"\t" + SpectrumParameter.getMaxSpecMzRange(charge));
			
			out.println("#ION\t" + charge + "\t" + iterationNum);
			for(IonType ion : sigIons){
				if(ion instanceof IonType.SuffixIon){
					out.println("s/"+ion.getCharge()+"/"+ion.getOffset() + "\t" + sigIonMap.get(ion));
				}else if(ion instanceof IonType.PrefixIon){
					out.println("p/"+ion.getCharge()+"/"+ion.getOffset() + "\t" + sigIonMap.get(ion));
				}else if(ion instanceof IonType.PrecursorIon){
					out.println("r/"+ion.getCharge()+"/"+ion.getOffset() + "\t" + sigIonMap.get(ion));
				}
			}
					
			out.println("#FEATURES");
			
			ArrayList<Feature> features = new ArrayList<Feature>();
			
			for(SpectrumParameter spar : independentFeatureMap.keySet()){
				HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> si = independentFeatureMap.get(spar);
				for(PeakParameter ppar : si.keySet()){
					HashMap<Integer, ArrayList<Feature>> ssi = si.get(ppar);
					for(int i=0; i<100; i++){
						ArrayList<Feature> cons = ssi.get(i);
						if(cons != null && !cons.isEmpty()){
				
							Collections.sort(cons, Collections.reverseOrder());
							int cntr = 0;
							for(Feature con:cons){
								if(!(con instanceof IntensityFeature) && cntr++ > minFeatureNumPerGroup && con.getDivergence() < KLDthreshold) break;
								features.add(con);
									
								out.println(con.toFileString());
								for(int j=0; j<sigIons.size(); j++){
																	out.print(con.getProbability(sigIons.get(j))+"\t");
								}								
								out.print(con.getProbability(IonType.NOISE)+"\t");
								out.println();
							
							}
						}
					}
				}
			}
			out.println("#############################################################");
			out.close();
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
}
