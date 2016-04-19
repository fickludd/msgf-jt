package uninovoOld.train;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.IonType;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;
import uninovoOld.IPDGenerator;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;

public class UniNovoTrainer {
	private String inputmgf;
	private Tolerance tol;
	private Tolerance pmtol;
	private int minCaseNum = 100;
	private int minSpecNum = 3000;
	private int maxConditionNum = 8000;
	private int specMzNum = 5;
	private int peakGroupNum = 10;
	private int peakRatioNum = 10;
	private int peakPartitionNum = 4;
	//private int peakPartition2Num = 5;
	
	private AminoAcidSet aaSet;
	private boolean discard = true;
	private int maxIonNum = 8;
	private boolean considerIonIntensity = false;
	private WindowFilter filter = null;
	
	public UniNovoTrainer(String inputmgf, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet){
		this.inputmgf = inputmgf;
		this.tol = tol;
		this.pmtol = pmtol;
		this.aaSet = aaSet;
	}
	
	public UniNovoTrainer filter(WindowFilter b) {filter = b; return this;}
	
	public UniNovoTrainer doNotDiscard() {discard = false; return this;}
	
	public UniNovoTrainer setMaxIonNum(int n) {maxIonNum = n;  return this;}
	
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
				//System.out.println(spar.getSpecMzRange() + "\t" + s.getPeptideMass());
				
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
	
	private float[] getMzRangeHistogram(int charge, int peakPartition2Num){
		Iterator<Spectrum> iterator;
		float[] numPerBin = new float[100];
		float sum = 0f;
		
		try {
			iterator = new SpectraIterator(inputmgf, new MgfSpectrumParser());
			while(iterator.hasNext()){
				Spectrum s = iterator.next();
				if(s.getCharge()!=charge) continue;
				numPerBin[Math.round(s.getPeptideMass()/121.6f)] ++;
				sum ++;
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
		sum /= peakPartition2Num;
		System.out.println("Mz range charge : " + charge);
	
		float[] partitionMzs = new float[peakPartition2Num-1];
		float k=0; int j=0;
		for(int i=0; i<numPerBin.length;i++){
			if(numPerBin[i]>0)
				System.out.println(i + "\t" +  numPerBin[i]);
			k+= numPerBin[i];
			if(k>=sum && j < partitionMzs.length){
				
				partitionMzs[j++]=i*121.6f;
				//System.out.println(i+ "\t" +k);
				k=0;
			}
		}

		
		return partitionMzs;
	}
	
	public void train(String para, int charge, int maxIterationNum, Enzyme enzyme){
		
		SpectrumParameter.setSpecMzRangeNum(specMzNum);
		PeakParameter.setGroupNum(peakGroupNum);
		PeakParameter.setPartitionNum(peakPartitionNum);
		PeakParameter.setPeakIntensityRatioNum(peakRatioNum);
		
		new File(para).delete();
		
		ArrayList<Integer> range = getChargeRange();
		
		for(int c=range.get(0); c<=Math.min(range.get(1), charge);c++){
			long time = System.currentTimeMillis();
			//SpectrumParameter.setPartitionMzs(getMzRangeHistogram(c, specMzNum), c);
			
			String mgf = new String(inputmgf);
			
			int maxMZR = getMaxSpecMzRange(c);
			SpectrumParameter.setMaxSpecMzRange(c, maxMZR);
			System.out.println("Charge : " + c + " Max Mz Range " + maxMZR);
			int sn =0;
			IonSelector ionfinder = new IonSelector(mgf, tol, filter, maxIonNum);
			sn = ionfinder.train(c, considerIonIntensity);
		
			for(int iterationNum = 0; iterationNum < maxIterationNum; iterationNum++){
				
				//if(iterationNum > 0){
				///	ionfinder.setSpecfilename(mgf);
				//	ionfinder.trainIonIonCurretRatio(c);
				//}
				

				ArrayList<IonType> sigIons = new ArrayList<IonType>();
				for(IonType ion : ionfinder.getSigIons().keySet()){
					if(!ion.equals(IonType.NOISE)) 
						sigIons.add(ion);
				}

				if(sigIons.isEmpty()) continue;//
				
				FeatureSelector confinder = new FeatureSelector(mgf, sigIons, tol, pmtol).filter(filter);
				confinder.train(c, iterationNum, aaSet);
				
				FeatureStatisticsTrainer contrainer = new FeatureStatisticsTrainer(mgf, ionfinder.getSigIons(), confinder.getSignificantConditions(), tol, pmtol).filter(filter);
				contrainer.train(c, iterationNum);
				contrainer.discardConditionsWithLessThan(minCaseNum, maxConditionNum, discard);
				contrainer.writeInFile(para, c, iterationNum);
				
				confinder = null;
				contrainer = null;
				
				
				String nextinputmgf = mgf.substring(0, mgf.length()-4) + iterationNum + ".mgf";
				
				if(iterationNum < maxIterationNum - 1){
					IPDGenerator gen = new IPDGenerator(para, aaSet, 0).filter(filter);
					gen.run(mgf, nextinputmgf, c, iterationNum+1, iterationNum, true, -1);
				}
				
				if(iterationNum>0) new File(mgf).delete();
				mgf =  nextinputmgf;
				IPDGenerator.clear();
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
	
	static private void printUsageAndExit(){
		System.out.println("Usage : java -cp UniNovo.jar uninovo.train.UniNovoTrainer" +
				"\n -i training mgf,mzXML, or ms2 file" +
				"\n -o output parameter file" +
				"\n -t ion tolerance (ending in ppm or Da)" +
				"\n -pt precursor ion tolerance (ending in ppm or Da)" 
				//"\n -aa amino acid file (if not specified, no modification except C+57 as a fixed mod will be considered.)"
				);
			
		System.out.println("Example : " +
				"\n java -Xmx2000m -cp UniNovo.jar uninovo.train.UniNovoTrainer -i /home/annotated.mgf -o /home/trained.par -pt 20ppm -t 50ppm");
		
		System.exit(0);
	}
	
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
			else if(arg.startsWith("-aa")) mode = 4;
			else{
				if(mode == 0) trainmgf = arg;
				else if(mode == 1) para = arg;
				else if(mode == 2) tol = Tolerance.parseToleranceStr(arg);
				else if(mode == 3) pmtol = Tolerance.parseToleranceStr(arg);
				else if(mode == 4) aaSet = AminoAcidSet.getAminoAcidSetFromModFile(arg);
			}
		}
		
		int maxIterationNum = 5;
		
		UniNovoTrainer trainer = new UniNovoTrainer(trainmgf, tol, pmtol, aaSet).filter(filter);//.setMaxIonNum(ionnum);
		trainer.train(para, maxCharge, maxIterationNum, enzyme);

	}
}
