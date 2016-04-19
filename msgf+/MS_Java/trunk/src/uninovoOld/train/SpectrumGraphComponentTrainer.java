package uninovoOld.train;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.IonType;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;
import parser.MgfSpectrumParser;
import uninovoOld.Edge;
import uninovoOld.IPDGenerator;
import uninovoOld.Node;
import uninovoOld.SpectrumGraph;
import uninovoOld.SpectrumGraphComponent;

public class SpectrumGraphComponentTrainer {

	private int maxAANum = 5;
	
	private Enzyme enzyme;
	private String specfilename, para;
	private Tolerance tol, pmtol;
	private AminoAcidSet aaSet;
	private WindowFilter filter;
	
	private float[] ionOffset;
	private float[][] ionWeights;
	
	//private float ionProbProbabilities[][][][];// casenum ion1 ion2 cur/false
	//private float enzymaticIonProbProbabilities[][][][];// casenum ion1 ion2 cur/false
	private float [][][][] edgeAccuracies; // minAANum, accL, accR, mass dev
	private float[] nodeNullAccuracies;
	private float[] nodeNoPeakAccuracies;
	private float[][] nodePosteriorProbabilities; // c/ic, minAANum
	private float[][][] edgeDeviationAccuracies; // c/ic, minAA, index
	/*
	private void trainNodeLR(int specCharge){
		Iterator<Spectrum> iterator;
		int length = 10;
		
		int sn = 0;
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPG ipe = new IPG(para, tol, pmtol, aaSet, 0).filter(filter);
			ArrayList<IonType> sigIons = ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(specCharge);

			ionProbProbabilities = new float[1<<sigIons.size()][length][length][2];// casenum ion1 ion2 cur/false
			enzymaticIonProbProbabilities = new float[1<<sigIons.size()][length][length][2];// casenum ion1 ion2 cur/false
			
			float[] ionProbProbabilitiesDivider = new float[2];
			      
			float[] enzymaticIonProbProbabilitiesDivider = new float[2];
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
			//	if(sn > 5000) break;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training node LR for NewSpectrumGraph : " + sn);
				}
				
				spec.setRanksOfPeaks();
		
				SpectrumGraph mag = new SpectrumGraph(spec, ipe, enzyme, tol, pmtol, false, true);
				
				ArrayList<Node> nodes = mag.getNodes();
				
				ArrayList<Float> prms = new ArrayList<Float>();
				
				for(float prm : spec.getAnnotation().getPRMMasses(true, 0)){
					prms.add(prm);
				}
			//	System.out.println(nodes.size());
				for(int i=1; i<nodes.size()-1;i++){
					Node node = nodes.get(i);
					int caseNumber = node.getCaseNumber();
					int bin1 = Math.round((length-1)*node.getIonProbsSelected()[0]); 
					int bin2 = Math.round((length-1)*node.getIonProbsSelected()[1]); 
				//	if(bin2 == 0) bin2 = bin1;
					
					boolean isCorrect = false;
					for(float cm : prms){
						if(node.isWithinTolerance(cm, spec.getPeptideMass())){
							isCorrect = true;
							break;
						}
					}
					
					if(node.isEnzymatic()){
						if(isCorrect){
							enzymaticIonProbProbabilities[caseNumber][bin1][bin2][0]++;
							enzymaticIonProbProbabilitiesDivider[0]++;
						}else{
							enzymaticIonProbProbabilities[caseNumber][bin1][bin2][1]++;
							enzymaticIonProbProbabilitiesDivider[1]++;
						}
					}else{
						if(isCorrect){
							ionProbProbabilities[caseNumber][bin1][bin2][0]++;
							ionProbProbabilitiesDivider[0]++;
						}else{
							ionProbProbabilities[caseNumber][bin1][bin2][1]++;
							ionProbProbabilitiesDivider[1]++;
						}
					}
				}
			}
			
			for(int c=0; c<enzymaticIonProbProbabilities.length; c++){
				for(int b1=0; b1<enzymaticIonProbProbabilities[c].length;b1++){
					for(int b2=0; b2<enzymaticIonProbProbabilities[c][b1].length;b2++){
						enzymaticIonProbProbabilities[c][b1][b2][0] /= enzymaticIonProbProbabilitiesDivider[0];
						enzymaticIonProbProbabilities[c][b1][b2][1] /= enzymaticIonProbProbabilitiesDivider[1];
					}
				}
			}
			
			for(int c=0; c<ionProbProbabilities.length; c++){
				for(int b1=0; b1<ionProbProbabilities[c].length;b1++){
					for(int b2=0; b2<ionProbProbabilities[c][b1].length;b2++){
						ionProbProbabilities[c][b1][b2][0] /= ionProbProbabilitiesDivider[0];
						ionProbProbabilities[c][b1][b2][1] /= ionProbProbabilitiesDivider[1];
					}
				}
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}*/
	
	
	private void trainNodeAccuracy(int specCharge){
		Iterator<Spectrum> iterator;
		double[][][] meanY;
		double[][][] meanYY;
		double[][][] meanXY;
		double[] meanX;
		double[] divider;
		float[] nullDivider, noPeakDivider;
		int sn = 0;
		
		SpectrumGraphComponent.setSectionNumber(SpectrumGraphComponent.maxSectionNum);
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPDGenerator ipe = new IPDGenerator(para, aaSet, 0).filter(filter);
			ArrayList<IonType> sigIons = ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(specCharge);
			
			nodePosteriorProbabilities = new float[2][maxAANum+1];
			
			nodeNullAccuracies = new float[SpectrumGraphComponent.maxSectionNum];
			nodeNoPeakAccuracies = new float[SpectrumGraphComponent.maxSectionNum];
			nullDivider = new float[SpectrumGraphComponent.maxSectionNum];
			noPeakDivider = new float[SpectrumGraphComponent.maxSectionNum];
			
			meanY = new double[1<<(sigIons.size())][][];
			meanYY = new double[1<<(sigIons.size())][][];
			meanXY = new double[1<<(sigIons.size())][][];
			meanX = new double[1<<(sigIons.size())];
			divider = new double[1<<(sigIons.size())];
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;

				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training node accuracy for SpectrumGraph : " + sn);
				}
				
				spec.setRanksOfPeaks();
		
				SpectrumGraph mag = new SpectrumGraph(spec, ipe, enzyme, tol, pmtol, false, true);
				
				ArrayList<Node> nodes = mag.getNodes();
				
				
				//nodeNoPeakAccuracy
				Random r = new Random();
				int num = 0;
				while(num++<200){
					float nm = r.nextFloat() * spec.getPeptideMass();
					Node nnode = new Node(nm, spec.getCharge(), null, tol, spec.getPeptideMass(), false, 0);

					boolean hasPeak = !mag.getNodesWithinTolerance(nm).isEmpty();
					if(!hasPeak){
						noPeakDivider[nnode.getSection()]++;
						//System.out.println(nnode.getSection() + "\t" + noPeakDivider[nnode.getSection()]);
					}
					
					if(nnode.isCorrect(spec.getAnnotation())){
						nodeNullAccuracies[nnode.getSection()]++;
						if(!hasPeak) nodeNoPeakAccuracies[nnode.getSection()]++;
					}
					
				
					nullDivider[nnode.getSection()]++;
				}
				
				HashSet<Node> cNodes = mag.getCorrectNodes(spec.getAnnotation());
				
				for(int i=1; i<nodes.size()-1;i++){
					Node node = nodes.get(i);
					int caseIndex=0;
					int cardinality = 0;
					float[] ionprobs = node.getIonProbs();
		
					int d=0;

					if(cNodes.contains(node)){
						d=1; 
					}
					
					if(node.getMinAA() < Integer.MAX_VALUE){
						if(cNodes.contains(node))
							nodePosteriorProbabilities[0][Math.min(nodePosteriorProbabilities[0].length-1, node.getMinAA())]++;
						else
							nodePosteriorProbabilities[1][Math.min(nodePosteriorProbabilities[1].length-1, node.getMinAA())]++;
						}
					
					for(int j=0; j<sigIons.size(); j++){ // ion
						if(ionprobs[j] > 0){
							caseIndex += 1<<j;
							cardinality++;
							if(cardinality >= Node.maxCardinality) break;
						}
					}
					
					if(cardinality==0){
						continue;
					}
					
					int[] indices = new int[cardinality];
					int t = 0;
					for(int j=0; j<sigIons.size(); j++){ // ion
						if(ionprobs[j] > 0){
							indices[t++] = j;
							if(t >= indices.length) break;
						}
					}
				
					divider[caseIndex]++;
					meanX[caseIndex] += d;
					
					if(meanY[caseIndex] == null){
						meanY[caseIndex] = new double[1][cardinality];
						meanYY[caseIndex] = new double[cardinality][cardinality];
						meanXY[caseIndex]  = new double[1][cardinality];
					}
					
					for(int j =0; j< cardinality;j++){ // ion
						float prob = ionprobs[indices[j]];
						meanY[caseIndex][0][j] += 	prob;
						meanXY[caseIndex][0][j] += d * prob;
						for(int k=0; k< cardinality; k++){
							float prob2 = ionprobs[indices[k]];
							meanYY[caseIndex][k][j] += prob2 * prob;
						}
					}
				}
			}
			
			for(int i=0; i<nodePosteriorProbabilities.length;i++){
				double sum = 0;
				for(int j=0; j<nodePosteriorProbabilities[i].length; j++){
					sum += nodePosteriorProbabilities[i][j];
				}
				for(int j=0; j<nodePosteriorProbabilities[i].length; j++){
					nodePosteriorProbabilities[i][j] /= sum;
				}
			}
			
			ionOffset= new float[1<<(sigIons.size())];
			ionWeights = new float[1<<(sigIons.size())][];

			for(int caseIndex = 1; caseIndex < 1<<(sigIons.size()); caseIndex++){
				if(divider[caseIndex] < 100) continue;
				
				meanX[caseIndex]/=divider[caseIndex];
			
				for(int i=0; i<meanY[caseIndex][0].length; i++){
					meanY[caseIndex][0][i]/=divider[caseIndex];
				}
				for(int i=0; i<meanXY[caseIndex][0].length; i++) meanXY[caseIndex][0][i] = 
					(meanXY[caseIndex][0][i]-divider[caseIndex]*meanX[caseIndex]*meanY[caseIndex][0][i])/(divider[caseIndex]-1);
				for(int i=0; i<meanYY[caseIndex].length; i++){
					for(int j=0; j<meanYY[caseIndex][i].length; j++) meanYY[caseIndex][i][j] = (meanYY[caseIndex][i][j] - divider[caseIndex]*meanY[caseIndex][0][i]*meanY[caseIndex][0][j])/(divider[caseIndex]-1);
				}
				
				meanYY[caseIndex] = MC.sum(meanYY[caseIndex], MC.multiply(MC.transpose(meanY[caseIndex]), meanY[caseIndex]), -1);
				meanYY[caseIndex] = MC.invert(meanYY[caseIndex]);
	
				double[][] coeff  = MC.multiply(meanYY[caseIndex], MC.transpose(MC.sum(meanXY[caseIndex], meanY[caseIndex], -meanX[caseIndex])));

				ionOffset[caseIndex] = (float) (meanX[caseIndex]-MC.multiply(meanY[caseIndex], coeff)[0][0]);
				ionWeights[caseIndex] = new float[coeff.length];
				for(int i=0; i<ionWeights[caseIndex].length;i++){
					ionWeights[caseIndex][i] = (float) coeff[i][0];
				}
			}
			
			for(int i=0; i<nullDivider.length; i++){
				//System.out.println(nodeNullAccuracies[i]  + "\t" + nullDivider[i]);
				nodeNullAccuracies[i] /= nullDivider[i];
				nodeNoPeakAccuracies[i] /= noPeakDivider[i];
			//	System.out.println(nodeNullAccuracies[i]);
			}
	}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
			
	
	/*
	private void trainIonProbConf(int specCharge, int maxRank, int maxIterationNum){
		Iterator<Spectrum> iterator;
		
		float[][][] divider; // sigIon * sigIon * 10 * 10
		int length = 10;
		
		int sn = 0;
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPG ipe = new IPG(para, tol, pmtol, aaSet).filter(filter);
			ArrayList<IonType> sigIons = ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(specCharge);

			ionProbCoeffs = new float[2][sigIons.size()][sigIons.size()][length][length];
			divider = new  float[2][sigIons.size()][sigIons.size()];
			
			while(iterator.hasNext()){
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				if(sn > 5000) break;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training ion prob coeffs for NewSpectrumGraph : " + sn);
				}
				
				spec.setRanksOfPeaks();
		
				ArrayList<NewNode> nodes = NewSpectrumGraph.getNodes(spec, ipe, para, enzyme, maxRank, tol, false);
				
				ArrayList<Float> prms = new ArrayList<Float>();
				for(float prm : spec.getAnnotation().getPRMMasses(true, 0)){
					prms.add(prm);
				}
				
				for(NewNode node : nodes){
					for(float cm : prms){
						if(node.isWithinTolerance(cm, spec.getPeptideMass(), tol)){
							node.isCorrect = true;
							break;
						}
					}
				}
				
				for(int i=1; i<nodes.size()-1;i++){
					NewNode node = nodes.get(i);
					int cur = node.isCorrect? 0 : 1;
					float[] ionprobs = node.getIonProbs();
		
					for(int j=0;j<ionprobs.length;j++){
						if(ionprobs[j] == 0) continue;
						//			int i2 = Math.round((divider[0].length-1)*rnode.getAccuracy()); 
						int bin1 = Math.round(ionprobs[j] * (length-1));
						for(int l=j;l<ionprobs.length;l++){
							if(ionprobs[l] == 0) continue;
							divider[cur][j][l]++;
							int bin2 = Math.round(ionprobs[l] * (length-1));
							ionProbCoeffs[cur][j][l][bin1][bin2]++;
						}
					}
				}
			}
			
			for(int k=0; k< ionProbCoeffs.length; k++){
				for(int j=0; j< ionProbCoeffs[k].length; j++){
					for(int l=j; l< ionProbCoeffs[k][j].length; l++){
						for(int x=0; x< ionProbCoeffs[k][j][l].length; x++){
							for(int y=0; y< ionProbCoeffs[k][j][l][x].length; y++){
								ionProbCoeffs[k][j][l][x][y] /= divider[k][j][l];
							}
						}
					}
				}
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
	private void trainEnzymaticNodeAccuracy(int specCharge, int maxRank, int maxIterationNum){
		if(enzyme == null) return;
		Iterator<Spectrum> iterator;
		
		int length = 20; // interpolation needed
		
		int sn = 0;
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPG ipg = new IPG(para, tol, pmtol, aaSet).filter(filter);
			accuracyForEnzymaticNodes = new float[length];
			float[] divider = new float[length];
			
			while(iterator.hasNext()){ 
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				if(sn > 5000) break;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training enyzmatic node accuracy : " + sn);
				}
		
			//	NewSpectrumGraph mag = new  NewSpectrumGraph(spec, ipg, para, maxRank, enzyme, tol, pmtol);
				
				ArrayList<NewNode> nodes = NewSpectrumGraph.getNodes(spec, ipg, para, enzyme, maxRank, pmtol, false);

				for(float prm : spec.getAnnotation().getPRMMasses(true, 0)){
					for(int i=1; i<nodes.size()-1;i++){
						NewNode node = nodes.get(i);
						if(node.isWithinTolerance(prm, spec.getPeptideMass(), tol)) node.isCorrect = true;
					}
				}
				
				for(int i=1; i<nodes.size()-1;i++){ // exclude source sink
					NewNode enode = nodes.get(i);
					if(enode.isEnzymaticNode()){
						if(enode.isNoPeakNode()) continue;
						int bin  = Math.round((accuracyForEnzymaticNodes.length-1)*enode.getAccuracy());
						if(enode.isCorrect)
							accuracyForEnzymaticNodes[bin]++;
						divider[bin]++;
					}
				}
			}
		
			for(int i=0;i<accuracyForEnzymaticNodes.length;i++){
				accuracyForEnzymaticNodes[i]/=divider[i];
			}
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
		ionWeights = null;
	}
	*/
	private void trainEdgeAccuracy(int specCharge){
		Iterator<Spectrum> iterator;
		
		int length = 20; // interpolation needed
		
		int sn = 0;
		
		SpectrumGraphComponent.setSectionNumber(SpectrumGraphComponent.maxSectionNum);
		
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			IPDGenerator ipg = new IPDGenerator(para, aaSet, 0).filter(filter);
			
			maxAANum = 5;
			float[][][][] divider = new float[maxAANum+1][length+1][length+1][1];
		//	float divider2 = 0;
			edgeAccuracies = new float[maxAANum+1][length+1][length+1][1];
			
			if(!SpectrumGraphComponent.useBinning(tol))
				edgeDeviationAccuracies = new float[2][maxAANum+1][Edge.massDeviationIndexNum];
		//	edgeNullAccuracies = 0;
			
			while(iterator.hasNext()){ 
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				//if(sn>1000)break;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training edge accuracy : " + sn);
				
				}
				
				SpectrumGraph mag = new SpectrumGraph(spec, ipg, enzyme, tol, pmtol, false, true);
			
				ArrayList<Node> nodes = mag.getNodes();
				HashSet<Edge> cEdges = mag.getCorrectEdges(spec.getAnnotation());
				
				for(int i=1; i<nodes.size()-1;i++){ // exclude sink
					Node rnode = nodes.get(i);
					
					/*Random r = new Random();
					int num = 0;
					while(rnode.isCorrect(spec.getAnnotation()) && num < 10){
						float m = rnode.getMass() + r.nextFloat() * (spec.getPeptideMass() - rnode.getMass());
						Node node = new Node(m, spec.getCharge(), null, ipg.getSigIonsOrderedByIntensityWithOutNoiseIon(specCharge), tol, false, 0);
						Edge edge = new Edge(rnode, node, spec.getPeptideMass(), pmtol, aaSet);
					//	if(edge.isValid()){
						
						if(edge.isCorrect(spec.getAnnotation())) edgeNullAccuracies ++;
						
						divider2 ++;
						num++;
					//	}
					} */
					//System.out.println(mag.getEdges(rnode));
					
					//nodePosteriorProbabilities acc, minAA
					int i2 = Math.round((divider[0].length-1)*rnode.getAccuracy()); 
					
					
					for(Edge e : mag.getLinkedEdges(rnode)){
						Node lnode = e.getLeftNode();
						
						if(lnode.getMass() == 0) continue;
					
						if(e.isPTM()) continue;
					//	if(e.getMinAANum() == Integer.MAX_VALUE) continue;
						int i1 = Math.min(divider.length-1, e.getMinAANum() );
						if(e.isCompositeMass()) i1 = 0;
						else if(i1 == 0) continue;
					
						
						int i3 = Math.round((divider[0][0].length-1)*lnode.getAccuracy()); 
							
						divider[i1][i2][i3][0] ++;
						
						if(cEdges.contains(e)){
							edgeAccuracies[i1][i2][i3][0] ++;					
						}
						
						if(edgeDeviationAccuracies!=null && e.getMassDeviationIndex()>=0) {
							int m = Math.min(edgeDeviationAccuracies[0].length-1, e.getMinAANum() );
							if(cEdges.contains(e)){
								edgeDeviationAccuracies[0][m][e.getMassDeviationIndex()]++;	
							}
							else{
								edgeDeviationAccuracies[1][m][e.getMassDeviationIndex()]++;	
							}
						}
					}
				}
			}
			
			for(int i=0;i<edgeAccuracies.length;i++){
				for(int j=0;j<edgeAccuracies[i].length;j++){
					for(int k=0;k<edgeAccuracies[i][j].length;k++){
						for(int l=0;l<edgeAccuracies[i][j][k].length;l++){
							if(divider[i][j][k][l] > 30){
								edgeAccuracies[i][j][k][l] /= divider[i][j][k][l];
							}else{
								edgeAccuracies[i][j][k][l] =0;
							}
								
						}
					}
				}
			}
			
	
		//	if(divider2 > 0){
		//		edgeNullAccuracies /= divider2;
		//	}
			
			if(edgeDeviationAccuracies!=null){
				for(int i=0;i<2;i++){
					for(int k=0; k< edgeDeviationAccuracies[i].length;k++){
						double sum = 0;
						for(int j=0; j<edgeDeviationAccuracies[i][k].length; j++){
							sum += edgeDeviationAccuracies[i][k][j];
						}
						for(int j=0; j<edgeDeviationAccuracies[i][k].length; j++){
							edgeDeviationAccuracies[i][k][j] /= sum;
							edgeDeviationAccuracies[i][k][j] = Math.max(1e-5f, edgeDeviationAccuracies[i][k][j]);
						}
					}
					
				}
				
				for(int i=0;i<2;i++){
					for(int k=0; k< edgeDeviationAccuracies[i].length;k++){
						double sum = 0;
						for(int j=0; j<edgeDeviationAccuracies[i][k].length; j++){
							sum += edgeDeviationAccuracies[i][k][j];
						}
						for(int j=0; j<edgeDeviationAccuracies[i][k].length; j++){
							edgeDeviationAccuracies[i][k][j] /= sum;
						}
					}
					
				}
			}
	
		}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
	
	private void write(String para, int charge){
		PrintWriter out;
		try {
			out = new PrintWriter(new FileWriter(para, true));
			
			out.println("##NODE&EDGEPARAMETERS##");
			if(ionWeights != null){
				int written = 0;

				out.println("#IONWEIGHT\t" + charge);
				for(int i=0; i< ionWeights.length;i++){
					boolean towrite = false;
			
					if(ionWeights[i] == null) continue;
					for(int j=0; j< ionWeights[i].length; j++){
						if(ionWeights[i][j] != 0) towrite = true;
					}
					if(!towrite) continue;
					if(ionOffset[i] > 0.3 || ionOffset[i] < -0.3) continue;
					written++;
					
					out.println("#OFF\t"+i + "\t" + ionOffset.length);
					out.println(ionOffset[i]);
					
					out.println("#WEIGHTS\t"+ i + "\t" +ionWeights[i].length+"\t" + ionWeights.length);
					for(float w : ionWeights[i]) out.println(w);
				}
				
			
				out.println("#ENDIONWEIGHT\t" + written);
				
				
			}
			
			if(nodeNullAccuracies != null){
				out.print("#NULLPROB\t" + charge);
				for(float p : nodeNullAccuracies)
					out.print("\t"+p);
				out.println();
				//out.println("#NOPEAKPROB\t" + charge + "\t" + nodeNoPeakAccuracy);
			}
			
			if(nodeNoPeakAccuracies != null){				
				out.print("#NOPEAKPROB\t" + charge);
				for(float p : nodeNoPeakAccuracies)
					out.print("\t"+p);
				out.println();
			}
			
			if(edgeDeviationAccuracies!=null){
				out.print("#DEVIATIONPROB1\t" + charge+"\t");
				for(int k=1; k<edgeDeviationAccuracies[0].length;k++){
					for(float n : edgeDeviationAccuracies[0][k]){
						out.print(n+":");
					}
					out.print("\t");
				}
				
				out.println();
				out.print("#DEVIATIONPROB2\t" + charge+"\t");
				for(int k=1; k<edgeDeviationAccuracies[1].length;k++){
					for(float n : edgeDeviationAccuracies[1][k]){
						out.print(n+":");
					}
					out.print("\t");
				}
				out.println();
			}
			
			if(nodePosteriorProbabilities != null){
				out.print("#NODEPOSTERIOR\t" + charge + "\t");
				for(int k=0;k<nodePosteriorProbabilities.length;k++){
					for(int j=1; j<nodePosteriorProbabilities[k].length; j++){
						out.print(nodePosteriorProbabilities[k][j]+":");
					}
					out.print("\t");
				}
				out.println();
			}
			
			/*if(accuracyForEnzymaticNodes!= null){
				out.println("#ENZYMATICNODESACCURACY\t" + charge+"\t" + accuracyForEnzymaticNodes.length);
				for(float acc : accuracyForEnzymaticNodes){
					out.println(acc);
				}
				out.println("#ENDENZYMATICNODESACCURACY\t" + charge);
			}*/
		
			if(edgeAccuracies != null){
				
				out.println("#EDGEACCURACY\t"+edgeAccuracies.length+"\t"+edgeAccuracies[0].length+"\t"+edgeAccuracies[0][0].length+"\t" + edgeAccuracies[0][0][0].length + "\t" +  charge);
				
				for(int i=0;i<edgeAccuracies.length;i++){
					out.println("##NUM\t"+i);
					for(int j=0;j<edgeAccuracies[i].length;j++){
						for(int k=0;k<edgeAccuracies[i][j].length;k++){
							for(int l=0;l<edgeAccuracies[i][j][k].length;l++){
								out.print(edgeAccuracies[i][j][k][l]+":");
							}
							out.print(" ");
							//System.out.println(edgeAccuracies[i][j][k] + "\t" + charge);
						}
						out.println();
					}
				
				}
				
				out.println("#ENDEDGEACCURACY");
			}
			
			/*if(edgeNullAccuracies != 0){
				
				out.println("#NULLEDGEACCURACY\t"+  charge);
				

				out.println("##NUM\t"+0);
				
				out.print(edgeNullAccuracies);
			
				out.println();
				
				out.println("#ENDNULLEDGEACCURACY");
			}*/
/*
			if(ionProbCoeffs != null){
				for(int c=0; c<ionProbCoeffs.length;c++){
					out.println("#IONPROBCOEFF\t"+c+"\t"+ionProbCoeffs[c].length+"\t"+ionProbCoeffs[c][0].length+"\t"+ionProbCoeffs[c][0][0].length+"\t" + ionProbCoeffs[c][0][0][0].length +"\t"+  charge);
					
					for(int i=0;i<ionProbCoeffs[c].length;i++){
						out.println("##NUM\t"+i);
						for(int j=0;j<ionProbCoeffs[c][i].length;j++){
							for(int k=0;k<ionProbCoeffs[c][i][j].length;k++){
								for(int l=0;l<ionProbCoeffs[c][i][j][k].length;l++){
									out.print(ionProbCoeffs[c][i][j][k][l]+":");
								}
								out.print(" ");
							}
							out.println();
						}
					
					}
					
					out.println("#ENDIONPROBCOEFF");
				}
			}
			*/
			/*if(ionProbProbabilities != null){
				out.println("#IONPROBPROB\t"+ionProbProbabilities.length+"\t"+ionProbProbabilities[0].length+"\t"+ionProbProbabilities[0][0].length+"\t" + ionProbProbabilities[0][0][0].length +"\t"+  charge);
				
				for(int i=0;i<ionProbProbabilities.length;i++){
					boolean toWrite = false;
					for(int j=0;j<ionProbProbabilities[i].length;j++){
						for(int k=0;k<ionProbProbabilities[i][j].length;k++){
							for(int l=0;l<ionProbProbabilities[i][j][k].length;l++){
								if(ionProbProbabilities[i][j][k][l]>0){
									toWrite = true;
								}
							}
						}
					}
					if(toWrite){
						out.println("##NUM\t"+i);
						for(int j=0;j<ionProbProbabilities[i].length;j++){
							for(int k=0;k<ionProbProbabilities[i][j].length;k++){
								for(int l=0;l<ionProbProbabilities[i][j][k].length;l++){
									out.print(ionProbProbabilities[i][j][k][l]+":");
								}
								out.print(" ");
							}
							out.println();
						}
					}
				}
				
				out.println("#ENDIONPROBPROB");
				
				out.println("#ENZIONPROBPROB\t"+enzymaticIonProbProbabilities.length+"\t"+enzymaticIonProbProbabilities[0].length+"\t"+enzymaticIonProbProbabilities[0][0].length+"\t" + enzymaticIonProbProbabilities[0][0][0].length +"\t"+  charge);
				
				for(int i=0;i<enzymaticIonProbProbabilities.length;i++){
					boolean toWrite = false;
					for(int j=0;j<enzymaticIonProbProbabilities[i].length;j++){
						for(int k=0;k<enzymaticIonProbProbabilities[i][j].length;k++){
							for(int l=0;l<enzymaticIonProbProbabilities[i][j][k].length;l++){
								if(enzymaticIonProbProbabilities[i][j][k][l]>0){
									toWrite = true;
								}
							}
						}
					}
					if(toWrite){
						out.println("##NUM\t"+i);
						for(int j=0;j<enzymaticIonProbProbabilities[i].length;j++){
							for(int k=0;k<enzymaticIonProbProbabilities[i][j].length;k++){
								for(int l=0;l<enzymaticIonProbProbabilities[i][j][k].length;l++){
									out.print(enzymaticIonProbProbabilities[i][j][k][l]+":");
								}
								out.print(" ");
							}
							out.println();
						}
					}
				}
				
				out.println("#ENDENZIONPROBPROB");
		
			}*/
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	static protected void checkParaFile(String para){
		try {
			File tempFile = new File(para + ".tmp");
			BufferedReader br;
			
			br = new BufferedReader(new FileReader(para));
		
			PrintWriter pw = new PrintWriter(new FileWriter(tempFile));
		
		    String s;
		    boolean readyToQuit =false;
		    while((s=br.readLine())!=null){
			    	if(readyToQuit){
		    		if(s.startsWith("##NODE&EDGEPARAMETERS##") || s.startsWith("#IONWEIGHT"))
		    			break;
		    	}
		    	
		    	pw.println(s);
		    	if(s.startsWith("#############################################################")){
		    		readyToQuit = true;
		    	}else
		    		readyToQuit = false;
		    }
		    
		    br.close();
		    pw.close();
		 //   System.exit(0);
		    new File(para).delete();
		    tempFile.renameTo(new File(para));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	protected void train(int charge, Enzyme enzyme){
		ionWeights = null;
		edgeAccuracies = null;
	//	ionProbProbabilities =null;
	
		this.enzyme = enzyme;
	//	IPG.setMaxIterationNum(0);
	//	trainNodeLR(charge, maxRank);
	//	write(para, charge);
		
	//	IPG.setMaxIterationNum(0);
		trainNodeAccuracy(charge);
		write(para, charge);
		
	//	trainNodeAccuracy(charge);
	//	write(para, charge);
	//	trainEnzymaticNodeAccuracy(charge, maxRank, maxIterationNum);
	//	write(para, charge);
		
	//	IPG.setMaxIterationNum(0);
		trainEdgeAccuracy(charge);
		write(para, charge);
	
	}
	
	protected SpectrumGraphComponentTrainer(String specfilename, String para, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet){
		this.specfilename = specfilename;
		this.para = para;
		this.tol = tol;
		this.pmtol = pmtol;
		this.aaSet = aaSet;
	}
	
	protected SpectrumGraphComponentTrainer filter(WindowFilter b) {filter = b; return this;}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String trainmgf = System.getProperty("user.home") + "/Dropbox/shewLengthAll.mgf"; //CID_Tryp_Confident.mgf Zubarev_HCD_Annotated.mgf shewLengthAll.mgf
		
		
		Enzyme enzyme = null;//Enzyme.TRYPSIN;
		
		//Tolerance tol = new Tolerance(10f, true);
		//Tolerance pmtol = new Tolerance(10f, true);
		Tolerance tol = new Tolerance(0.5f, false);
		Tolerance pmtol =  new Tolerance(20f, true);
		WindowFilter filter = new WindowFilter(6, 50);
		
		if(args.length > 0) trainmgf = args[0];
		if(args.length > 1) {
			if(args[1].contains("h")){
				tol = new Tolerance(20f, true);
			}
		}
		if(args.length > 2) {
			if(args[2].contains("L")){
				//enzyme = Enzyme.LysC;
			}
		}
		
		
		String para = System.getProperty("user.home") +"/Dropbox/PAR" + trainmgf.substring(trainmgf.lastIndexOf("/"), trainmgf.lastIndexOf(".")) + ".pars";
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		UniNovoTrainer trainer = new UniNovoTrainer(trainmgf, tol, pmtol, aaSet).filter(filter);
		
		SpectrumGraphComponentTrainer.checkParaFile(para);
		ArrayList<Integer> range = trainer.getChargeRange();
		
		for(int c=range.get(0); c<=Math.min(range.get(1), 100);c++){
			SpectrumGraphComponentTrainer mtrainer = new SpectrumGraphComponentTrainer(trainmgf, para, tol, pmtol,  aaSet).filter(filter);
			mtrainer.train(c, enzyme);
			mtrainer = null;
		}

	}

}


