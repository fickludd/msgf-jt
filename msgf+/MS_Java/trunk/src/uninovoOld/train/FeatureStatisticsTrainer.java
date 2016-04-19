package uninovoOld.train;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.IonType;
import msutil.Peak;
import msutil.Spectrum;
import msutil.WindowFilter;
import uninovoOld.IPDGenerator;
import uninovoOld.UniNovo;
import uninovoOld.features.Feature;
import uninovoOld.features.IntensityFeature;
import uninovoOld.features.LinkingFeature;
import uninovoOld.features.OffsetFeature;
import uninovoOld.parameters.PeakParameter;
import uninovoOld.parameters.SpectrumParameter;

public class FeatureStatisticsTrainer {
	static int minFeatureNumPerGroup = 50;// tmp TODO
	static float KLDthreshold = 0;
	
	private String specfilename;
	private Tolerance tol;
	private Tolerance pmtol;
	private ArrayList<IonType> sigIons;
	private HashMap<IonType, Float> sigIonMap;
	private HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> conditionMap;
	private HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>> independentConditionMap;
	private HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>> rankIonProbMap;
	
	private WindowFilter filter = null;
	
	public FeatureStatisticsTrainer(String specfilename,  HashMap<IonType, Float> sigIonMap, HashMap<SpectrumParameter, HashMap<PeakParameter, ArrayList<Feature>>> conditionMap,
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
		
		this.conditionMap = conditionMap;
		this.tol = tol;
		this.pmtol = pmtol;
		this.rankIonProbMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>>();
	}
	
	public FeatureStatisticsTrainer filter(WindowFilter b) {filter = b; return this;}
	
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
			int maxRank = IPDGenerator.getMaxRank(spec);
			
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
			HashMap<PeakParameter, ArrayList<Feature>> sigconMap = conditionMap.get(spar);
			if(!rankIonProbMap.containsKey(spar))
				rankIonProbMap.put(spar, new HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>>());
			
			HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>> srankIonProbMap = rankIonProbMap.get(spar);
			
			for(Peak bp : spec){
				if(bp.getRank() > maxRank) continue;
				
				PeakParameter ppar = new PeakParameter(bp, spec,iterationNum);
				ArrayList<Feature> sigcons = sigconMap.get(ppar);
				ArrayList<Feature> matchedcons = new ArrayList<Feature>();
				IntensityFeature nullCondition = null;
				IonType expIon = null;
				
				if(!srankIonProbMap.containsKey(ppar))
					srankIonProbMap.put(ppar, new HashMap<Integer, HashMap<IonType, Float>>());
				HashMap<Integer, HashMap<IonType, Float>> ssr = srankIonProbMap.get(ppar);
				
				if(!ssr.containsKey(bp.getRank()))
					ssr.put(bp.getRank(), new HashMap<IonType, Float>());
				
				HashMap<IonType, Float> sssr = ssr.get(bp.getRank());
				
				for(IonType ion : sigIons){
					if(pgen.isExplainedBy(bp, ion, tol, pmtol)){
						expIon = ion;
						break;
					}
				}
				
				if(expIon == null) expIon = IonType.NOISE;
				
				if(!sssr.containsKey(expIon))
					sssr.put(expIon, 0f);
				
				sssr.put(expIon, sssr.get(expIon)+1);
				
				for(Feature con : sigcons){
					if(con.holdsFor(bp, filteredspec, spar, ppar, tol, pmtol, iterationNum)){
						matchedcons.add(con);
					}
					if(con instanceof IntensityFeature){
						nullCondition = (IntensityFeature) con;			
					}	
				}
				
				for(Feature con : matchedcons){
					con.registerNullCondition(nullCondition);
					con.addIonCount(expIon);
				}
			}
		}
	}
	
	protected void discardConditionsWithLessThan(int minCaseNum, int maxConditionNum){
		discardConditionsWithLessThan(minCaseNum, maxConditionNum, true);
	}
	
	public void discardConditionsWithLessThan(int minCaseNum, int maxConditionNum, boolean discard){
		ArrayList<Feature> toErase = new ArrayList<Feature>();
		
		int[] num = new int[4];
		int[] erasednum = new int [3];
		ArrayList<Float> KLDs = new ArrayList<Float>();
		
		for(SpectrumParameter spar : conditionMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			for(PeakParameter ppar : cm.keySet()){
				for(Feature con : cm.get(ppar)){
					if(con instanceof IntensityFeature){
						num[0]++;
						con.calculateIonProbs();
					}
				}
			}
		}
	
		for(SpectrumParameter spar : conditionMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			for(PeakParameter ppar : cm.keySet()){
				ArrayList<Feature> cs = cm.get(ppar);
				for(int i=0; i<cs.size();i++){
					Feature con = cs.get(i);
					if(con instanceof IntensityFeature) continue;
					else if(con instanceof OffsetFeature) num[1]++;
					else if(con instanceof LinkingFeature) num[2]++;
			//		else if(con instanceof DensityFeature) num[3]++;
					
					float sum = con.calculateIonProbs();
					if(discard && sum < minCaseNum){
						cs.remove(i--);
						toErase.add(con);
					}else{
						KLDs.add(con.getKLDivergenceFromIntensityFeature());
					}
				}
			}
		}
		
		Collections.sort(KLDs, Collections.reverseOrder());
		KLDthreshold = 0;

		if(!KLDs.isEmpty())
			KLDthreshold = KLDs.get(Math.min(KLDs.size()-1, maxConditionNum));
		
		for(SpectrumParameter spar : conditionMap.keySet()){
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			for(PeakParameter ppar : cm.keySet()){
				ArrayList<Feature> cs = cm.get(ppar);
				for(int i=0; i<cs.size();i++){
					Feature con = cs.get(i);
					if(con instanceof IntensityFeature) continue;
				
					if(discard && con.getKLDivergenceFromIntensityFeature() < KLDthreshold){
					//	toErase.add(con);
					//	cs.remove(i--);
					}
			
				}
			}
		}
		
		
		for(Feature con : toErase){
			if(con instanceof OffsetFeature) erasednum[0]++;
			else if(con instanceof LinkingFeature) erasednum[1]++;
	//		else if(con instanceof DensityFeature) erasednum[2]++;
		}
		
		System.out.println("# Null feature : " + num[0]);
		System.out.println("# Ion dependency feature : " + num[1] + " -> " + (num[1] - erasednum[0]));
		System.out.println("# Linking feature : " + num[2] + " -> " + (num[2] - erasednum[1]));
	//	System.out.println("# Density Condition : " + num[3] + " -> " + (num[3] - erasednum[2]));
		
	}
	
	private void updateIndependentConditionMap(){
		/* 0 = null
		 * 1 = bridging
		 * 2 = density
		 * 3 = gof normal
		 * 4 ~ = gof with charge off, comp, diff base peak charges*/
		
		independentConditionMap = new HashMap<SpectrumParameter, HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>>();

		for(SpectrumParameter spar : conditionMap.keySet()){
			int sc = spar.getSpecCharge();
			
			HashMap<PeakParameter, ArrayList<Feature>> cm = conditionMap.get(spar);
			
			if(!independentConditionMap.containsKey(spar)) independentConditionMap.put(spar,  new HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>>());
			HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> si = independentConditionMap.get(spar);
				
			for(PeakParameter ppar : cm.keySet()){
				if(!si.containsKey(ppar)) si.put(ppar, new HashMap<Integer, ArrayList<Feature>>());
				HashMap<Integer, ArrayList<Feature>> ssi = si.get(ppar);
				for(Feature con : cm.get(ppar)){
					int id;
					if(con instanceof IntensityFeature) id = 0;
					else if(con instanceof LinkingFeature){
						id = con.getBasePeakCharge();
						//id = 1;
					}else{
						OffsetFeature gcon = (OffsetFeature) con;
						/*if(gcon.getChargeOffset() == 0){
							if(!gcon.isComplementary()) id = 2;
							else id = 3;
						}else{
							id = 4 + sc * (gcon.getBasePeakCharge() - 1) + (gcon.getBasePeakCharge() + gcon.getChargeOffset() - 1);
						}*/
						
						if(!gcon.isComplementary() && gcon.getChargeOffset() == 0) id = sc + 1;
						else{
							id = sc + 2 + (gcon.isComplementary()? sc * sc : 0) + sc * (gcon.getBasePeakCharge() - 1) + (gcon.getBasePeakCharge() + gcon.getChargeOffset() - 1);
						}
					}
					if(!ssi.containsKey(id)) ssi.put(id, new ArrayList<Feature>());
					ArrayList<Feature> cons = ssi.get(id);
					
					cons.add(con);
				}
			}	
		}
	}
	
	public void writeInFile(String file, int charge, int iterationNum){
		writeInFile(file, charge, iterationNum, false);
	}
	
	public void writeInFile(String file, int charge, int iterationNum, boolean readable){
		updateIndependentConditionMap();
		
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
			
			//if(sigIonMap.containsKey(IonType.NOISE));
			//	out.println("n/"+0+"/"+0 + "\t" + sigIonMap.get(IonType.NOISE));
			//out.println("#END");
		//	}
			
			/*out.println("#RANKPROB\t"+iterationNum+"\t" + charge);
			
			int max = 0, min = Integer.MAX_VALUE;
			for(SpectrumParameter spar : rankIonProbMap.keySet()){
				HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>> sr = rankIonProbMap.get(spar);
				for(PeakParameter ppar : sr.keySet()){
					HashMap<Integer, HashMap<IonType, Float>> ssr = sr.get(ppar);
					for(int i : ssr.keySet()){
						max = Math.max(max,i);
						min = Math.min(min, i);
						
						float sum = 0;
						for(IonType ion : ssr.get(i).keySet()){
							sum += ssr.get(i).get(ion);
						}
						for(IonType ion : ssr.get(i).keySet()){
							ssr.get(i).put(ion, ssr.get(i).get(ion)/sum);
						}
					}
				}
			}
		
			for(SpectrumParameter spar : rankIonProbMap.keySet()){
				out.println("#SPAR\t"+spar.toFileString());
				HashMap<PeakParameter, HashMap<Integer, HashMap<IonType, Float>>> sr = rankIonProbMap.get(spar);
				for(PeakParameter ppar : sr.keySet()){
					out.println("#PPAR\t"+ppar.toFileString());
					HashMap<Integer, HashMap<IonType, Float>> ssr = sr.get(ppar);
					for(int i=min; i<= max; i++){
						if(!ssr.containsKey(i)) continue;
						out.print(i+"\t");
						for(IonType ion : sigIons){
							Float p = ssr.get(i).get(ion);
							if(p == null) p = 0f;
							out.print(p+"\t");
						}
						out.println();
					}
				}
			}*/
			out.println("#FEATURES");
			
			ArrayList<Feature> features = new ArrayList<Feature>();
		//	FeatureComparator comp = new FeatureComparator(tol);
			
			for(SpectrumParameter spar : independentConditionMap.keySet()){
				HashMap<PeakParameter, HashMap<Integer, ArrayList<Feature>>> si = independentConditionMap.get(spar);
				for(PeakParameter ppar : si.keySet()){
					HashMap<Integer, ArrayList<Feature>> ssi = si.get(ppar);
					for(int i=0; i<100; i++){
						ArrayList<Feature> cons = ssi.get(i);
						if(cons != null && !cons.isEmpty()){
							if(!readable)
								out.println("#"+i);
							
							Collections.sort(cons, Collections.reverseOrder());
							int cntr = 0;
							for(Feature con:cons){
								if(!(con instanceof IntensityFeature) && cntr++ > minFeatureNumPerGroup && con.getKLDivergenceFromIntensityFeature() < KLDthreshold) break;
							//	if(cntr++ > maxFeatureNum) break;
								
								
								features.add(con);
									
								if(!readable){
									out.println(con.toFileString());// + "$" + con.getKLDivergenceFromNullCondition() + "$$"+con.getNullCondition());
									for(int j=0; j<sigIons.size(); j++){
										if(con instanceof IntensityFeature)
											out.print(con.getIonProbability(sigIons.get(j))+"\t");
										else
											out.print(con.getProbability(sigIons.get(j))+"\t");
									}
									if(con instanceof IntensityFeature)
										out.print(con.getIonProbability(IonType.NOISE)+"\t");
									else
										out.print(con.getProbability(IonType.NOISE)+"\t");
									out.println();
								}
							}
						}
					}
				}
			}
			if(readable){
				
				Collections.sort(features, Collections.reverseOrder());
				
				
				
				
//				for(Feature f : features){
//					if(f instanceof IntensityFeature){
//						IntensityFeature ff = (IntensityFeature)f;
//						int t = 10-ff.getBasePeakParameter().getBasePeakGroupNum();
//						
//						if(t!=10 || f.getBasePeakCharge()!=1) continue;
//						out.print("0" + "&" + "&"+"&"+"&"+"1"+"&"+"&");
//						float sum = 1;
//	
//						for(int i=0;i<sigIons.size();i++){
//							float prob = ff.getIntensityFeature().getIonProbability(sigIons.get(i));
//							sum -= prob;
//							out.print(String.format("%.0f&", prob*100f));
//						}
//						out.println(String.format("%.0f&0.0", sum*100f));
//					}
//				}
				
				int cntr = 1;
				out.print("#&t&x&r&T&z_1&z_2");
				for(int i=0;i<sigIons.size();i++){
					out.print("&");
				}
				out.println("");
				for(Feature f : features){
					if(f instanceof OffsetFeature){
						OffsetFeature ff = (OffsetFeature)f;
						if(!ff.isPresent()) continue;
						int t = 10-ff.getBasePeakParameter().getBasePeakGroupNum();
						
						//if(t!=10) continue;
						float x = ff.getIOFF().getOffset();
						int r = ff.getPeakIntensityRatio();
						int c1 = ff.getBasePeakCharge();
						int c2 = c1 + ff.getIOFF().getChargeOffset();
						int T=0;
						
						if(ff.getIOFF().isComplementary() || c1 != c2) continue;
						float sum = 1;
						
						if(ff.getIOFF().isComplementary()) T = 1;
						out.print(cntr++ + "&" +t+ "&" + String.format("%.2f", x) + "&"+r+"&"+T+"&"+c1+"&"+c2+"&   ");
						
						for(int i=0;i<sigIons.size();i++){
							float prob = ff.getIonProbability(sigIons.get(i));
							sum -= prob;
							out.print(String.format("%.0f&", prob*100f));
						}
						out.println(String.format("%.0f", sum*100f));
						if(cntr>10) break;
						
					}
				
				}
			
				cntr = 1;
				out.print("#&t&x&r&T&z_1&z_2");
				for(int i=0;i<sigIons.size();i++){
					out.print("&");
				}
				out.println("");
				for(Feature f : features){
					if(f instanceof OffsetFeature){
						OffsetFeature ff = (OffsetFeature)f;
						if(!ff.isPresent()) continue;
						int t = 10-ff.getBasePeakParameter().getBasePeakGroupNum();
						
						//if(t!=10) continue;
						float x = ff.getIOFF().getOffset();
						int r = ff.getPeakIntensityRatio();
						int c1 = ff.getBasePeakCharge();
						int c2 = c1 + ff.getIOFF().getChargeOffset();
						int T=0;
						if(!ff.getIOFF().isComplementary() || c1 != c2) continue;
						float sum = 1;
				
						if(ff.getIOFF().isComplementary()) T = 1;
						out.print(cntr++ + "&" +t+ "&" + String.format("%.2f", x) + "&"+r+"&"+T+"&"+c1+"&"+c2+"&   ");
						
						for(int i=0;i<sigIons.size();i++){
							float prob = ff.getIonProbability(sigIons.get(i));
							sum -= prob;
							out.print(String.format("%.0f&", prob*100f));
						}
						out.println(String.format("%.0f", sum*100f));
						if(cntr>10) break;
					}
				
				}
			
				cntr = 1;
				out.print("#&t&x&r&T&z_1&z_2");
				for(int i=0;i<sigIons.size();i++){
					out.print("&");
				}
				out.println("");
				for(Feature f : features){
					if(f instanceof OffsetFeature){
						OffsetFeature ff = (OffsetFeature)f;
						if(!ff.isPresent()) continue;
						int t = 10-ff.getBasePeakParameter().getBasePeakGroupNum();
						
						//if(t!=10) continue;
						float x = ff.getIOFF().getOffset();
						int r = ff.getPeakIntensityRatio();
						int c1 = ff.getBasePeakCharge();
						int c2 = c1 + ff.getIOFF().getChargeOffset();
						int T=0;
						
						if(ff.getIOFF().isComplementary() || c1 == c2) continue;
						float sum = 1;
				
						if(ff.getIOFF().isComplementary()) T = 1;
						out.print(cntr++ + "&" +t+ "&" + String.format("%.2f", x) + "&"+r+"&"+T+"&"+c1+"&"+c2+"&   ");
						
						for(int i=0;i<sigIons.size();i++){
							float prob = ff.getIonProbability(sigIons.get(i));
							sum -= prob;
							out.print(String.format("%.0f&", prob*100f));
						}
						out.println(String.format("%.0f", sum*100f));
						if(cntr>10) break;
					}
				
				}
			
				cntr = 1;
				out.print("#&t&x&r&T&z_1&z_2");
				for(int i=0;i<sigIons.size();i++){
					out.print("&");
				}
				out.println("");
				for(Feature f : features){
					if(f instanceof OffsetFeature){
						OffsetFeature ff = (OffsetFeature)f;
						if(!ff.isPresent()) continue;
						int t = 10-ff.getBasePeakParameter().getBasePeakGroupNum();
						
						//if(t!=10) continue;
						float x = ff.getIOFF().getOffset();
						int r = ff.getPeakIntensityRatio();
						int c1 = ff.getBasePeakCharge();
						int c2 = c1 + ff.getIOFF().getChargeOffset();
						int T=0;
						if(!ff.getIOFF().isComplementary() || c1 == c2) continue;
						float sum = 1;
				
						if(ff.getIOFF().isComplementary()) T = 1;
						out.print(cntr++ + "&" +t+ "&" + String.format("%.2f", x) + "&"+r+"&"+T+"&"+c1+"&"+c2+"&   ");
						
						for(int i=0;i<sigIons.size();i++){
							float prob = ff.getIonProbability(sigIons.get(i));
							sum -= prob;
							out.print(String.format("%.0f&", prob*100f));
						}
						out.println(String.format("%.0f", sum*100f));
						if(cntr>10) break;
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
