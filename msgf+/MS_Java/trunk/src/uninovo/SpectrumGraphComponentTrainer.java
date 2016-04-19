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

import uninovo.parser.MgfSpectrumParser;
import uninovo.util.AminoAcidSet;
import uninovo.util.Enzyme;
import uninovo.util.IonType;
import uninovo.util.SpectraIterator;
import uninovo.util.Spectrum;
import uninovo.util.Tolerance;
import uninovo.util.WindowFilter;

/**
 * The Class SpectrumGraphComponentTrainer: a part of UniNovoTrainer that learns edge and node statistics
 */
public class SpectrumGraphComponentTrainer {

	/**
	 * Check if the parameter file is the last version (in terms of iteration number).
	 *
	 * @param para the parameter file name
	 */
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
		    new File(para).delete();
		    tempFile.renameTo(new File(para));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/** The aa set. */
	private AminoAcidSet aaSet;
	
	/** The edge accuracies. */
	private float [][][][] edgeAccuracies; // minAANum, accL, accR, mass dev
	
	/** The edge deviation accuracies. */
	private float[][][] edgeDeviationAccuracies; // c/ic, minAA, index
	
	/** The enzyme. */
	private Enzyme enzyme;
	
	/** The filter. */
	private WindowFilter filter;
	
	/** The ion offset. */
	private float[] ionOffset;
	
	/** The ion weights. */
	private float[][] ionWeights;
	
	/** The max aa num. */
	private int maxAANum = 5;
	
	/** The node no peak accuracies. */
	private float[] nodeNoPeakAccuracies;
	
	/** The node null accuracies. */
	private float[] nodeNullAccuracies;
	
	/** The node posterior probabilities. */
	private float[][] nodePosteriorProbabilities; // c/ic, minAANum
	
	/** The para. */
	private String specfilename, para;
	
	/** The pmtol. */
	private Tolerance tol, pmtol;
	
	/**
	 * Instantiates a new spectrum graph component trainer.
	 *
	 * @param specfilename the spec file name
	 * @param para the parameter file name
	 * @param tol the MS2 tol
	 * @param pmtol the MS1 tol
	 * @param aaSet the amino acid set
	 */
	protected SpectrumGraphComponentTrainer(String specfilename, String para, Tolerance tol, Tolerance pmtol, AminoAcidSet aaSet){
		this.specfilename = specfilename;
		this.para = para;
		this.tol = tol;
		this.pmtol = pmtol;
		this.aaSet = aaSet;
	}
	
	/**
	 * Filter.
	 *
	 * @param b the window filter
	 * @return the spectrum graph component trainer
	 */
	protected SpectrumGraphComponentTrainer filter(WindowFilter b) {filter = b; return this;}
	
	/**
	 * Train.
	 *
	 * @param charge the charge
	 * @param enzyme the enzyme
	 */
	protected void train(int charge, Enzyme enzyme){
		ionWeights = null;
		edgeAccuracies = null;	
		this.enzyme = enzyme;

		trainNodeAccuracy(charge);
		write(para, charge);
	
		trainEdgeAccuracy(charge);
		write(para, charge);
	
	}
	
	/**
	 * Train edge accuracy.
	 *
	 * @param specCharge the charge
	 */
	private void trainEdgeAccuracy(int specCharge){
		Iterator<Spectrum> iterator;
		
		int length = 20; // interpolation needed
		
		int sn = 0;
		
		SpectrumGraphComponent.setMassGroupNumber(SpectrumGraphComponent.maxMassGroupNum);
		
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			PeakIonProbabilityGenerator ipg = new PeakIonProbabilityGenerator(para, aaSet, 0).setFilter(filter);
			
			maxAANum = 5;
			float[][][][] divider = new float[maxAANum+1][length+1][length+1][1];
			edgeAccuracies = new float[maxAANum+1][length+1][length+1][1];
			
			if(!SpectrumGraphComponent.useBinning(tol))
				edgeDeviationAccuracies = new float[2][maxAANum+1][Edge.massDeviationIndexNum];
	
			while(iterator.hasNext()){ 
				Spectrum spec = iterator.next();
				
				if(specCharge > 0 && spec.getCharge() != specCharge) continue;
				if(spec.getAnnotation().isModified()) continue;
				
				sn++;
				if(prevsn < sn/1000){
					prevsn = sn/1000;
					System.out.println("training edge accuracy : " + sn);
				
				}
				
				SpectrumGraph mag = new SpectrumGraph(spec, ipg, enzyme, tol, pmtol, false, true);
			
				ArrayList<Node> nodes = mag.getFilteredNodes();
				HashSet<Edge> cEdges = mag.getCorrectEdges(spec.getAnnotation());
				
				for(int i=1; i<nodes.size()-1;i++){ // exclude sink
					Node rnode = nodes.get(i);
			
					int i2 = Math.round((divider[0].length-1)*rnode.getAccuracy()); 
					
					
					for(Edge e : mag.getLinkedEdges(rnode)){
						Node lnode = e.getLeftNode();
						
						if(lnode.getMass() == 0) continue;
					
						if(e.isPTM()) continue;
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
	
	/**
	 * Train node accuracy.
	 *
	 * @param specCharge the charge
	 */
	private void trainNodeAccuracy(int specCharge){
		Iterator<Spectrum> iterator;
		double[][][] meanY;
		double[][][] meanYY;
		double[][][] meanXY;
		double[] meanX;
		double[] divider;
		float[] nullDivider, noPeakDivider;
		int sn = 0;
		
		SpectrumGraphComponent.setMassGroupNumber(SpectrumGraphComponent.maxMassGroupNum);
		
		try {
			iterator = new SpectraIterator(specfilename, new MgfSpectrumParser());
			int prevsn = 0; 
			PeakIonProbabilityGenerator ipe = new PeakIonProbabilityGenerator(para, aaSet, 0).setFilter(filter);
			ArrayList<IonType> sigIons = ipe.getSigIonsOrderedByIntensityWithOutNoiseIon(specCharge);
			
			nodePosteriorProbabilities = new float[2][maxAANum+1];
			
			nodeNullAccuracies = new float[SpectrumGraphComponent.maxMassGroupNum];
			nodeNoPeakAccuracies = new float[SpectrumGraphComponent.maxMassGroupNum];
			nullDivider = new float[SpectrumGraphComponent.maxMassGroupNum];
			noPeakDivider = new float[SpectrumGraphComponent.maxMassGroupNum];
			
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
				
				ArrayList<Node> nodes = mag.getFilteredNodes();
				
				
				//nodeNoPeakAccuracy
				Random r = new Random();
				int num = 0;
				while(num++<200){
					float nm = r.nextFloat() * spec.getPeptideMass();
					Node nnode = new Node(nm, spec.getCharge(), null, tol, spec.getPeptideMass(), false, 0);

					boolean hasPeak = !mag.getNodesWithinTolerance(nm).isEmpty();
					if(!hasPeak){
						noPeakDivider[nnode.getLocationGroup()]++;
					}
					
					if(nnode.isCorrect(spec.getAnnotation())){
						nodeNullAccuracies[nnode.getLocationGroup()]++;
						if(!hasPeak) nodeNoPeakAccuracies[nnode.getLocationGroup()]++;
					}
					
				
					nullDivider[nnode.getLocationGroup()]++;
				}
				
				HashSet<Node> cNodes = mag.getCorrectNodes(spec.getAnnotation());
				
				for(int i=1; i<nodes.size()-1;i++){
					Node node = nodes.get(i);
					int caseIndex=0;
					int cardinality = 0;
					float[] ionprobs = node.getFPV();
		
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
							if(cardinality >= Node.maxIonNumberPerNode) break;
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
				nodeNullAccuracies[i] /= nullDivider[i];
				nodeNoPeakAccuracies[i] /= noPeakDivider[i];
			}
	}catch (FileNotFoundException e) {
			System.exit(0);
			e.printStackTrace();
		} 
	}
	
	/**
	 * Write the learned statistics on a file
	 *
	 * @param para the parameter file name
	 * @param charge the charge
	 */
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
						}
						out.println();
					}
				
				}
				
				out.println("#ENDEDGEACCURACY");
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}


