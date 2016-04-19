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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import uninovo.parser.MgfSpectrumParser;
import uninovo.util.AminoAcidSet;
import uninovo.util.Enzyme;
import uninovo.util.IonType;
import uninovo.util.SpectraIterator;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
import uninovo.util.WindowFilter;

/**
 * The Class UniNovoTrainer learns all statistics used in UniNovo
 */
public class UniNovoTrainer {
	
	/**
	 * The main method.
	 *
	 * @param args the arguments
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	static public void main(String[] args) throws IOException{

		if(args.length < 1) printUsageAndExit();
		
		int maxCharge = 100;
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol = new Tolerance(20f, true);
		WindowFilter filter = new WindowFilter(6, 50);
		Enzyme enzyme = null;

		String trainmgf = "";
		String para = "";
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		int mode = -1;
		for(String arg : args){
			if(arg.startsWith("-i")) mode = 0;
			else if(arg.startsWith("-o")) mode = 1;
			else if(arg.startsWith("-t")) mode = 2;
			else if(arg.startsWith("-pt")) mode = 3;
			else{
				if(mode == 0) trainmgf = arg;
				else if(mode == 1) para = arg;
				else if(mode == 2) tol = Tolerance.parseToleranceStr(arg);
				else if(mode == 3) pmtol = Tolerance.parseToleranceStr(arg);
			}
		}
		
		int maxIterationNum = 5;
		
		UniNovoTrainer trainer = new UniNovoTrainer(trainmgf, tol, pmtol, aaSet).filter(filter);//.setMaxIonNum(ionnum);
		trainer.train(para, maxCharge, maxIterationNum, enzyme);

	}
	
	/**
	 * Prints the usage and exit.
	 */
	static private void printUsageAndExit(){
		System.out.println("Usage : java -cp UniNovo.jar uninovo.train.UniNovoTrainer" +
				"\n -i training mgf,mzXML, or ms2 file" +
				"\n -o output parameter file" +
				"\n -t ion tolerance (ending in ppm or Da)" +
				"\n -pt precursor ion tolerance (ending in ppm or Da)" 
				);
			
		System.out.println("Example : " +
				"\n java -Xmx2000m -cp UniNovo.jar uninovo.train.UniNovoTrainer -i /home/annotated.mgf -o /home/trained.par -pt 20ppm -t 50ppm");
		
		System.exit(0);
	}
	
	/** The amino acid set. */
	private AminoAcidSet aaSet;
	
	/** If ion intensity is considered or not. If not, ion current values will be considered. Under test */
	private boolean considerIonIntensity = false;
	
	/** The discard - discard features with low divergences */
	private boolean discard = true;
	
	/** The window filter. */
	private WindowFilter filter = null;
	
	/** The input mgf. */
	private String inputmgf;
	
	/** The max feature num - total number of features per charge should not exceed this number */
	private int maxFeatureNum = 8000;
	
	/** The max ion num - max number of ion types */
	private int maxIonNum = 8;
	
	/** The min case num - features appeared less than this number will be discarded. */
	private int minCaseNum = 100;
	
	/** The min spec num - the number of spectra of a specific charge should exceed this number */
	private int minSpecNum = 3000;
	
	/** The peak group num. see PeakParameter*/
	private int peakGroupNum = 10;
	
	/** The peak partition num. see PeakParameter*/
	private int peakPartitionNum = 4;
	
	/** The peak ratio num. see PeakParameter*/
	private int peakRatioNum = 10;
	
	/** The MS1 tol. */
	private Tolerance pmtol;
	
	/** The number of groups of spectra according to their m/z values. */
	private int specMzNum = 5;
	
	/** The MS2 tol. */
	private Tolerance tol;
	
	/**
	 * Instantiates a new UniNovo trainer.
	 *
	 * @param inputmgf the input mgf
	 * @param tol the MS2 tol
	 * @param pmtol the MS1 tol
	 * @param aaSet the aa set
	 */
	public UniNovoTrainer(String inputmgf, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet){
		this.inputmgf = inputmgf;
		this.tol = tol;
		this.pmtol = pmtol;
		this.aaSet = aaSet;
	}
	
	/**
	 * Do not discard.
	 *
	 * @return this
	 */
	public UniNovoTrainer doNotDiscard() {discard = false; return this;}
	
	/**
	 * Filter.
	 *
	 * @param b the window filter
	 * @return this
	 */
	public UniNovoTrainer filter(WindowFilter b) {filter = b; return this;}
	
	
	/**
	 * Gets the charge range.
	 *
	 * @return the charge range
	 */
	ArrayList<Integer> getChargeRange(){
		Iterator<Spectrum> iterator;
		int[] numPerCharge = new int[100];
		int min = 100, max = 0;
		
		try {
			iterator = new SpectraIterator(inputmgf, new MgfSpectrumParser());
			while(iterator.hasNext()){
				Spectrum s = iterator.next();
				int charge = s.getCharge();
				
				numPerCharge[charge]++;

			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
		
		for(int i = 0 ; i<numPerCharge.length; i++){
			if(numPerCharge[i] > minSpecNum){
				max = max < i ? i  : max;
				min = min > i ? i : min;
			}
		}
		
		ArrayList<Integer> range = new ArrayList<Integer>();
		range.add(min); range.add(max);
		System.out.println("Charge range : " + min + "-" + max);
		return range;
		
		
	}
	
	/**
	 * Gets the max spec mz range.
	 *
	 * @param charge the charge
	 * @return the max spec mz range
	 */
	private int getMaxSpecMzRange(int charge){
		Iterator<Spectrum> iterator;
		int[] n = new int[100];
		try {
			iterator = new SpectraIterator(inputmgf, new MgfSpectrumParser());
			while(iterator.hasNext()){
				Spectrum s = iterator.next();
				if(s.getCharge()!=charge) continue;
				SpectrumParameter spar = new SpectrumParameter(s);
				n[spar.getSpecMzRange()]++;			
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
		
		int maxn = 0;
		for(int i=n.length-1; i>=0; i--){
			if(n[i]>150){
				maxn = i;
				break;
			}
		}
		return maxn;
	}
	
	/**
	 * Sets the max ion num.
	 *
	 * @param n the n
	 * @return this
	 */
	public UniNovoTrainer setMaxIonNum(int n) {maxIonNum = n;  return this;}
	
	/**
	 * Train.
	 *
	 * @param para the parameter file name
	 * @param charge the charge
	 * @param maxIterationNum the max iteration number
	 * @param enzyme the enzyme
	 */
	public void train(String para, int charge, int maxIterationNum, Enzyme enzyme){
		
		SpectrumParameter.setSpecMzRangeNum(specMzNum);
		PeakParameter.setGroupNum(peakGroupNum);
		PeakParameter.setPartitionNum(peakPartitionNum);
		PeakParameter.setPeakIntensityRatioNum(peakRatioNum);
		
		new File(para).delete();
		
		ArrayList<Integer> range = getChargeRange();
		
		for(int c=range.get(0); c<=Math.min(range.get(1), charge);c++){
			long time = System.currentTimeMillis();
			String mgf = new String(inputmgf);
			
			int maxMZR = getMaxSpecMzRange(c);
			SpectrumParameter.setMaxSpecMzRange(c, maxMZR);
			System.out.println("Charge : " + c + " Max Mz Range " + maxMZR);
			int sn =0;
			IonSelector ionfinder = new IonSelector(mgf, tol, filter, maxIonNum);
			sn = ionfinder.train(c, considerIonIntensity);
		
			for(int iterationNum = 0; iterationNum < maxIterationNum; iterationNum++){
				ArrayList<IonType> sigIons = new ArrayList<IonType>();
				for(IonType ion : ionfinder.getSigIons().keySet()){
					if(!ion.equals(IonType.NOISE)) 
						sigIons.add(ion);
				}

				if(sigIons.isEmpty()) continue;//
				
				FeatureSelector confinder = new FeatureSelector(mgf, sigIons, tol, pmtol).filter(filter);
				confinder.train(c, iterationNum, aaSet);
				
				FeatureStatisticsTrainer contrainer = new FeatureStatisticsTrainer(mgf, ionfinder.getSigIons(), confinder.getSignificantFeatures(), tol, pmtol).filter(filter);
				contrainer.train(c, iterationNum);
				contrainer.discardFeaturesWithLowDivergences(minCaseNum, maxFeatureNum, discard);
				contrainer.writeInFile(para, c, iterationNum);
				
				confinder = null;
				contrainer = null;
				String nextinputmgf = mgf.substring(0, mgf.length()-4) + iterationNum + ".mgf";
				
				if(iterationNum < maxIterationNum - 1){
					PeakIonProbabilityGenerator gen = new PeakIonProbabilityGenerator(para, aaSet, 0).setFilter(filter);
					gen.updateSpectrumFile(mgf, nextinputmgf, c, iterationNum);
				}
				
				if(iterationNum>0) new File(mgf).delete();
				mgf =  nextinputmgf;
				PeakIonProbabilityGenerator.clear();
			}
			System.out.println("Charge " + c + " training Done : " + (float)(System.currentTimeMillis() - time)/sn/1000 + " sec/spec");
			new File(mgf).delete();
			ionfinder = null;
		}
		
		
		SpectrumGraphComponentTrainer.checkParaFile(para);
		
		for(int c=range.get(0); c<=Math.min(range.get(1), charge);c++){
			SpectrumGraphComponentTrainer mtrainer = new SpectrumGraphComponentTrainer(inputmgf, para, tol, pmtol,  aaSet).filter(filter);
			mtrainer.train(c, enzyme);
			mtrainer = null;
		}
	}
}
