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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import uninovo.parser.BufferedLineReader;
import uninovo.parser.MgfSpectrumParser;
import uninovo.util.AminoAcidSet;
import uninovo.util.IonType;
import uninovo.util.Peak;
import uninovo.util.SpectraIterator;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
import uninovo.util.WindowFilter;


/**
 * Calculate the peak ion probability (i.e., the probability that a peak represents a specific ion peak)
 * that will be used to calculate FPVs for each prefix mass (or a node in a spectrum graph). In the paper, 
 * we described that FPVs are calculated directly from features; in practice, ion probabilities are easier to 
 * deal with. FPVs are then obtained from the ion probabilities in SpectrumGraph class.
 * @author kyowon
 */

public class PeakIonProbabilityGenerator {
	
	//hash table for independent feature groups
	//key : fragmentation method index, iteration number, spectrum parameter, peak parameter, group index 
	//value: features
	/** The independent feature map. */
	static  private HashMap<Integer, HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>> independentFeatureMap = null;
	
	//hash table for ions
	//key : fragmentation method index, charge, ion  
	//value: average intensity
	/** The ion charge map. */
	static private HashMap<Integer, HashMap<Integer, HashMap<IonType, Float>>>  ionChargeMap = null;
	
	//parent correction mode
	/** The is for pm correction. */
	static private boolean isForPMCorrection = false;
	
	//MS2 tolerance
	/** The pmtols. */
	static private Tolerance[] pmtols;
	
	//hash table for ions sorted according to the descending order of average intensities
	//key : fragmentation method index, charge
	//value: sorted ions
	/** The sorted ion charge map. */
	static private HashMap<Integer, HashMap<Integer, ArrayList<IonType>>> sortedIonChargeMap = null;
	
	//MS1 tolerances
	/** The tols. */
	static private Tolerance[] tols;
	
	/**
	 * Clears all static data - used only for training
	 */
	static public void clear(){
		ionChargeMap = null;
		sortedIonChargeMap = null; 
		independentFeatureMap = null;
	}
	
	/**
	 * Get maximum intensity rank to be used.
	 *
	 * @param spec spectrum
	 * @return the max rank
	 */
	static public int getMaxRank(Spectrum spec){ return Math.round(spec.getPeptideMass()/121.6f * 20);}//Integer.MAX_VALUE;}
	
	/**
	 * Check if the current mode is for parent mass correction
	 * @return the mode
	 */
	public static boolean isForPMCorrection() {
		return isForPMCorrection;
	}
	
	/**
	 * Change the mode of pips generation for parent mass correction
	 * @param isForPMCorrection set for parent mass correction mode
	 */
	public static void setForPMCorrection(boolean isForPMCorrection) {
		PeakIonProbabilityGenerator.isForPMCorrection = isForPMCorrection;
	}
	
	//amino acid set
	/** The aa set. */
	private AminoAcidSet aaSet;
	
	//window filter
	/** The filter. */
	private WindowFilter filter = null;
	
	//fragmentation method index
	/** The fragmentation method index. */
	private int fragmentationMethodIndex = 0;
	
	//max charge trained 
	/** The max charge. */
	private int maxCharge = 100;
	
	//max iteration number
	/** The max iteration number. */
	private int maxIterationNumber = 0;	

	//min charge trained
	/** The min charge. */
	private int minCharge = 2;
	
	//parameter file name
	/** The parafile. */
	private String parafile;
	
	/**
	 * Constructor.
	 *
	 * @param parafile the name of the parameter file
	 * @param aaSet the set of amino acids (required to define linking features)
	 * @param fragmentationMethodIndex the fragmentation method index
	 */
	public PeakIonProbabilityGenerator(String parafile, AminoAcidSet aaSet, int fragmentationMethodIndex){
		this.parafile = parafile;
		this.aaSet = aaSet;
		this.fragmentationMethodIndex = fragmentationMethodIndex;
		readFromFile();
		maxCharge = getBoundedChargeRange().get(1);
		minCharge = getBoundedChargeRange().get(0);
	}
		
	/**
	 * Returns the amino acid set
	 * @return amino acid set 
	 */
	public AminoAcidSet getAASet() {return aaSet;}
	
	/**
	 * Returns the range of bounded charges (i.e., charges trained)
	 * @return the range of bounded charges 
	 */
	private ArrayList<Integer> getBoundedChargeRange(){
		ArrayList<Integer> r = new ArrayList<Integer>();
		int min = 1000, max = 0;
		for(int c : ionChargeMap.get(fragmentationMethodIndex).keySet()){
			min = Math.min(min, c);
			max = Math.max(max, c);
		}
		r.add(min); r.add(max);
		return r;
	}
	
	/**
	 * Get the maximum charge trained
	 * @return the maximum charge
	 */
	public int getMaxCharge() {
		return maxCharge;
	}
	
	/**
	 * Get the maximum iteration number
	 * @return the maximum iteration number
	 */
	private int getMaxIterationNum() {
		return maxIterationNumber;
	}
	
	/**
	 * Get the minimum charge trained
	 * @return the minimum charge
	 */
	public int getMinCharge() {
		return minCharge;
	}
	
	/**
	 * Returns pips (that is obtained after the last iteration)
	 * @param spec spectrum
	 * @return pips for peaks in the spectrum after the last iteration
	 */
	public HashMap<Peak, HashMap<IonType, Float>> getPIPs(Spectrum spec){ 
		HashMap<Peak, HashMap<IonType, Float>> pips = null;
		 
		int maxIterationNum = Math.min(getMaxIterationNum(), readFromFile());
		
		if(isForPMCorrection()) 
			maxIterationNum = 1;		
		
		int maxRank = getMaxRank(spec);				
		for(int iterationNum = 0; iterationNum < maxIterationNum; iterationNum++){
			 pips = getPIPsForEachIteration(spec, maxRank, fragmentationMethodIndex, iterationNum);	
			 updateSpectrum(spec, pips, fragmentationMethodIndex);			 
		}
		return pips;	
	}
	
	/**
	 * Calculate pips for an iteration.
	 *
	 * @param spec spectrum
	 * @param maxRank that maximum rank of peak to be considered
	 * @param fragmentationMethodIndex the fragmentation method index
	 * @param iterationNumber iteration number
	 * @return pips for peaks in the spectrum
	 */
	private HashMap<Peak, HashMap<IonType, Float>> getPIPsForEachIteration(Spectrum spec, int maxRank, int fragmentationMethodIndex, int iterationNumber){
		HashMap<Peak, HashMap<IonType, Float>> pips = new HashMap<Peak, HashMap<IonType, Float>>();	
		spec.setRanksOfPeaks();		
		Spectrum filteredspec = spec;		
		if(filter != null) filteredspec = filter.apply(spec);						
		SpectrumParameter spar = new SpectrumParameter(spec);		
		HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> sigFeatureMap = independentFeatureMap.get(fragmentationMethodIndex).get(iterationNumber).get(spar);
		if(sigFeatureMap == null || sigFeatureMap.isEmpty()) return pips;
		
		PeakParameter.refresh();
		
		for(Peak bp : spec){
			if(bp.getRank() > maxRank) continue;
			PeakParameter ppar = new PeakParameter(bp, spec, iterationNumber);
			
			if(ppar.getBasePeakGroupNum() > PeakParameter.getMaxGroupNum()) continue;

			HashMap<Integer, ArrayList<Feature>> featureMap = sigFeatureMap.get(ppar);		
			ArrayList<Feature> satisfiedFeatures = new ArrayList<Feature>();
			
			for(int setnum : featureMap.keySet()){	
				ArrayList<Feature> features = featureMap.get(setnum);
				
				//assert(!features.isEmpty());
				for(int i = 0; i< features.size(); i++){
					Feature feature = features.get(i);
					if(feature instanceof LinkingFeature) continue;
					if(isForPMCorrection()){						
						if(feature instanceof OffsetFeature){
							OffsetFeature ofeature = (OffsetFeature)feature;
							if(!ofeature.isComplementary()) continue;
						}else if(!(feature instanceof IntensityFeature))
							continue;
					}
					if(feature.isSatisfiedBy(bp, filteredspec, spar, ppar, tols[fragmentationMethodIndex], pmtols[fragmentationMethodIndex], iterationNumber)){
						satisfiedFeatures.add(feature);						
						break;
					}			
				}
			}

			if(!satisfiedFeatures.isEmpty())
				pips.put(bp, Feature.getPeakIonProbabilityDistribution(bp, satisfiedFeatures, getSigIonsOrderedByIntensityWithOutNoiseIon(spar.getSpecCharge()), spec));
			
		}
		return pips;
	}
	
	/**
	 * Returns the range of bounded charges (i.e., charges trained)
	 *
	 * @param charge the charge
	 * @return the range of bounded charges
	 */
	public ArrayList<IonType> getSigIonsOrderedByIntensityWithOutNoiseIon(int charge) {
		charge = Math.min(charge, maxCharge);
		charge = Math.max(charge, minCharge);
		
		if(ionChargeMap.get(fragmentationMethodIndex) == null) return new ArrayList<IonType>();
		
		if(sortedIonChargeMap == null) sortedIonChargeMap = new HashMap<Integer, HashMap<Integer, ArrayList<IonType>>>();
		
		if(!sortedIonChargeMap.containsKey(fragmentationMethodIndex)){
			sortedIonChargeMap.put(fragmentationMethodIndex, new HashMap<Integer, ArrayList<IonType>>());
		}
		
		HashMap<Integer, ArrayList<IonType>>
			subsortedIonChargeMap = sortedIonChargeMap.get(fragmentationMethodIndex);

		if(!subsortedIonChargeMap.containsKey(charge))
			subsortedIonChargeMap.put(charge, new ArrayList<IonType>());
		
		
		ArrayList<IonType> ret = subsortedIonChargeMap.get(charge);		
		if(ret != null && !ret.isEmpty()) return ret;		
		ArrayList<Float> intensities = new ArrayList<Float>();	
		for(IonType ion : ionChargeMap.get(fragmentationMethodIndex).get(charge).keySet()){
			intensities.add(ionChargeMap.get(fragmentationMethodIndex).get(charge).get(ion));
		}
		
		Collections.sort(intensities);
		
		for(int i = intensities.size()-1; i>=0; i--){
			for(IonType ion : ionChargeMap.get(fragmentationMethodIndex).get(charge).keySet()){
				if(intensities.get(i) == ionChargeMap.get(fragmentationMethodIndex).get(charge).get(ion)){
					ret.add(ion);
				}
			}		
		}
		ret.remove(IonType.NOISE);		
		return ret;
	}

	/**
	 * Reads data from parameter file
	 * @return the maximum iteration number found in the parameter file
	 */
	private int readFromFile(){
		if(getMaxIterationNum() > 0) return getMaxIterationNum();

		BufferedLineReader in;
		int maxIterationNum = 0;
		System.out.println("Reading parameters for type index " + fragmentationMethodIndex + " : " +  parafile);

		try {
			in = new BufferedLineReader(parafile);
			String s;
			int mode = -1;
			int groupIndex = -1;
			int chargeForIon =  -1;
			HashMap<IonType, Float> ionProbMap = null;
			Feature feature = null;
			Feature nullFeature = null;
			
			HashMap<Integer, ArrayList<IonType>>
				sigIonMap = new HashMap<Integer, ArrayList<IonType>>();
			
			while((s = in.readLine()) != null){
			
				if(s.startsWith("#ION\t")){
					chargeForIon = Integer.parseInt(s.split("\t")[1]);
					if(ionChargeMap == null)
						ionChargeMap = new HashMap<Integer,HashMap<Integer, HashMap<IonType, Float>>>();
					if(!ionChargeMap.containsKey(fragmentationMethodIndex))
						ionChargeMap.put(fragmentationMethodIndex, new HashMap<Integer, HashMap<IonType, Float>>());

					mode = 0; continue;
				}
				if(s.equals("#FEATURES")){
					if(independentFeatureMap == null)
						independentFeatureMap = new HashMap<Integer, HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>>();
					if(!independentFeatureMap.containsKey(fragmentationMethodIndex))
						independentFeatureMap.put(fragmentationMethodIndex, new HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>());

					mode = 1; 
					continue;
				}
				if(s.startsWith("#SPECPARTITION\t")){
					chargeForIon = Integer.parseInt(s.split("\t")[1]);
					mode = 2;
					continue;
				}			
				
				if(s.startsWith("#PMTOL")){
					if(pmtols == null) pmtols = new Tolerance[100];
					pmtols[fragmentationMethodIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
					mode = -1;
					continue;
				}
				if(s.startsWith("#TOL")){
					if(tols == null) tols = new Tolerance[100];
					tols[fragmentationMethodIndex] = Tolerance.parseToleranceStr(s.split("\t")[1]);
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
						continue;
					}					
	
					HashMap<Integer, HashMap<IonType, Float>> ionMap = ionChargeMap.get(fragmentationMethodIndex);
					
					if(!ionMap.containsKey(chargeForIon)){
						ionMap.put(chargeForIon, new HashMap<IonType, Float>());
					}
					
					HashMap<IonType, Float> ions = ionMap.get(chargeForIon);
				
					if(!token[0].startsWith("n"))
						ions.put(IonType.getIonType(token[0]), Float.parseFloat(token[1]));
					else ions.put(IonType.NOISE, Float.parseFloat(token[1]));
			
					if(!sigIonMap.containsKey(chargeForIon))
						sigIonMap.put(chargeForIon, new ArrayList<IonType>());
					
					ArrayList<IonType> is = sigIonMap.get(chargeForIon);
					
					IonType iont = IonType.getIonType(token[0]);
					if(!is.contains(iont)) is.add(iont);

				}else if(mode == 1){
					if(s.startsWith("#")){
						groupIndex = Integer.parseInt(s.substring(1));
					}else if(s.startsWith("N")){
						feature = IntensityFeature.parseFileString(s);
						nullFeature = feature;
					}else if(s.startsWith("G")){
						feature = OffsetFeature.parseFileString(s);
					}else if(s.startsWith("B")){
						feature = LinkingFeature.parseFileString(s, aaSet);				
					}else{
						ionProbMap = new HashMap<IonType, Float>();
						String[] token = s.split("\t");
						
						ArrayList<IonType> ions = sigIonMap.get(feature.getSpectrumParameter().getSpecCharge());
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
						
						feature.setIonProbMap(ionProbMap);
						feature.setIntensityFeature(nullFeature);
						updateIndependentFeatureMap(feature, groupIndex);
						maxIterationNum = Math.max(feature.getIterationNum() + 1, maxIterationNum);
						
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
				}
			}
			
			in.close();			
			SpectrumGraphComponent.read(parafile, fragmentationMethodIndex);
		}catch (IOException e) {
			System.exit(1);
			e.printStackTrace();
		}
		setMaxIterationNum(maxIterationNum);
		
		System.out.println("Iteration Number : " + maxIterationNum);
		return maxIterationNum;
	}

	/**
	 * Set a window filter
	 * @param b filter
	 * @return this 
	 */
	public PeakIonProbabilityGenerator setFilter(WindowFilter b) {filter = b; return this;}

	/**
	 * Set the maximum iteration number
	 * @param maxIterationNum the maximum iteration number to be set
	 */
	private void setMaxIterationNum(int maxIterationNum) {
		this.maxIterationNumber = maxIterationNum;
	}
	
	/**
	 * Put a feature into one of independent feature groups
	 * @param feature feature
	 * @param groupIndex the index of group
	 */
	private void updateIndependentFeatureMap(Feature feature, int groupIndex){
		int iterationNum = feature.getIterationNum();
		SpectrumParameter spar = feature.getSpectrumParameter();
		PeakParameter ppar = feature.getBasePeakParameter();
		
		if(!independentFeatureMap.containsKey(fragmentationMethodIndex))
			independentFeatureMap.put(fragmentationMethodIndex, new HashMap<Integer, HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>>());
			
		if(!independentFeatureMap.get(fragmentationMethodIndex).containsKey(iterationNum))
			independentFeatureMap.get(fragmentationMethodIndex).put(iterationNum, new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>());
		HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>
			s1 = independentFeatureMap.get(fragmentationMethodIndex).get(iterationNum);
		
		if(!s1.containsKey(spar))
			s1.put(spar, new HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>());
		HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>
			s2 = s1.get(spar);
	
		if(!s2.containsKey(ppar))
			s2.put(ppar, new HashMap<Integer, ArrayList<Feature>>());
		HashMap<Integer, ArrayList<Feature>> s3 = s2.get(ppar);
	
		if(!s3.containsKey(groupIndex))
			s3.put(groupIndex, new ArrayList<Feature>());
		ArrayList<Feature> features = s3.get(groupIndex);
		//
		int index = Collections.binarySearch(features, feature, Collections.reverseOrder());
		if(index < 0) index = -index-1;
		features.add(index, feature);
	
	}

	/**
	 * Update the intensities of peaks in a spectrum according to their pips.
	 * @param spec spectrum
	 * @param pips pips
	 * @param fragmentationMethodIndex fragmentationMethodIndex
	 */
	private void updateSpectrum(Spectrum spec, HashMap<Peak, HashMap<IonType, Float>> pips, int fragmentationMethodIndex){
		 float high = 0, low = 0;		
		 SpectrumParameter spar = new SpectrumParameter(spec);
	     HashMap<IonType, Float> ionfactorMap = ionChargeMap.get(fragmentationMethodIndex).get(spar.getSpecCharge());
	
		 for(Peak p : spec){
			 if(p.getRank() == 1){
				 high = p.getIntensity();
				 break;
			 }
		}		
	
		 low = high*1e-2f;
		 
		 for(Peak p : pips.keySet()){
			HashMap<IonType, Float> probs = pips.get(p);
			float intensity = 0; 
				
			float factor = 1;
			
			ArrayList<IonType> ions = getSigIonsOrderedByIntensityWithOutNoiseIon(spec.getCharge());
			
			for(int i=0; i<ions.size(); i++){
				IonType ion = ions.get(i);
				if(ion instanceof IonType.PrecursorIon) continue;
				if(ion.equals(IonType.NOISE)) continue;						
				factor = ionfactorMap.get(ion);				
				intensity += probs.get(ion) *  factor;

			}
			intensity = intensity * (high - low) + low;
			p.setIntensity(intensity);
		 }		 		
	}

	/**
	 * Update the spectra according to their pips; pips are calculated inside. Only used for training
	 * @param specfilename spectrum file name
	 * @param outfilename output spectrum file name
	 * @param specCharge target charge for spectra
	 * @param iterationNumber iteration number that is trained
	 */	
	public void updateSpectrumFile(String specfilename, String outfilename, int specCharge, int iterationNumber){
		
		long time = System.currentTimeMillis();
		
		setMaxIterationNum(0); // init
		int maxIterationNum = Math.min(iterationNumber+1, readFromFile());
		setMaxIterationNum(maxIterationNum); // update
		
		int sn = 0;
		try {
			PrintStream out = null;
			if(outfilename != null){
				out = new PrintStream(outfilename);				
			}
			
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
				sn++;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("running " +": " + sn);
				}
										
				spec.setRanksOfPeaks();
				int maxRank = Math.round(spec.getPeptideMass()/1000*150);
				
				for(int l=0; l< spec.size();l++){ 
					Peak p = spec.get(l);
					if(p.getRank() >= maxRank){ 
						spec.remove(p);
						l--;
					}
				}
				HashMap<Peak, HashMap<IonType, Float>> pips = null;
				for(int iterationNum = iterationNumber; iterationNum < maxIterationNum; iterationNum++){
					 pips = getPIPsForEachIteration(spec, maxRank, fragmentationMethodIndex, iterationNum);
					 updateSpectrum(spec, pips, fragmentationMethodIndex);
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
					
					spec.outputMgf(out);
				}				
			}
			out.close();
		}catch (FileNotFoundException e) {
			System.exit(1);
			e.printStackTrace();
		}
		
		System.out.println("Running Done : " + (float)(System.currentTimeMillis() - time)/sn/1000 + " sec/spec");
		
	}
	
}
