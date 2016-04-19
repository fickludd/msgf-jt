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

import uninovo.util.AminoAcid;
import uninovo.util.Peptide;
import uninovo.util.Tolerance;

/**
 * The Class Node represents a node in a spectrum graph
 */
public class Node extends SpectrumGraphComponent{
	
	/** The Constant max ion number per node. */
	public final static int maxIonNumberPerNode = 3;
	
	/**
	 * calculate l rfrom accuracy.
	 *
	 * @param nacc the null accuracy - accuracy of node without supporting ion.
	 * @param acc the accuracy
	 * @param isNoPeak the is no peak
	 * @return the float
	 */
	static private float calculateLRfromAccuracy(double nacc, double acc){
		double up = 1-nacc;		
	
		acc = getBoundedValue(acc, up, nacc);
		return (float) (Math.log(acc/(1-acc) * (1-nacc)/nacc));
	}
	
	/**
	 * Calculate mass group.
	 *
	 * @param mass the mass
	 * @param pm the parent mass
	 * @return the mass group
	 */
	static private int calculateMassGroup(float mass, float pm){
		return Math.max(0, Math.min(massGroupNumber-1, (int)(mass/pm * massGroupNumber)));
	}
	
	/** The attachment cost used to calculate accuracy*/
	private float attachmentCost = Float.MAX_VALUE;
	
	/** The ion probs. */
	private float[] fpv;
	
	/** The is empty enzymatic node meaning an enzymatic node without any supporting ion */
	protected boolean isEmptyEnzymaticNode = false;
	
	/** is enzymatic node? */
	private boolean isEnzymaticNode = false;
	
	/** is sink. */
	private boolean isSink = false;
	
	/** The mass location group - rougly describes where is this node in terms of its mass*/
	private int locationGroup = 0;
	
	/** The min aa. */
	private int minAA = Integer.MAX_VALUE;
	
	/** The posterior accuracy - accuracy after considering edge accuracy. */
	private float posteriorAccuracy = -1;
	
	/** The posterior lr - lr after considering edge lr */
	private float posteriorLR=-Float.MAX_VALUE;
	
	/**
	 * Instantiates a new node.
	 *
	 * @param mass the mass
	 * @param charge the charge
	 * @param fpv the ion probs
	 * @param tol the MS2 tol
	 * @param pm the parent mass
	 * @param isEnzymaticNode is enzymatic node
	 * @param fragmentationMethodIndex the fragmentationMethodIndex
	 */
	public Node(float mass, int charge, float[] fpv, Tolerance tol, float pm, boolean isEnzymaticNode, int fragmentationMethodIndex){
		super(mass, charge, fragmentationMethodIndex);
		this.fpv = fpv;		
		this.isEnzymaticNode = isEnzymaticNode;
		this.locationGroup = calculateMassGroup(mass, pm);
		this.setTol(tol);
		setLRandAccuracy();
	}
	
	
	/**
	 * Instantiates a new node.
	 *
	 * @param mass the mass
	 * @param tol the MS2 tol
	 * @param pm the parent mass
	 * @param fragmentationMethodIndex the fragmentationMethodIndex
	 */
	public Node(float mass, Tolerance tol, float pm, int fragmentationMethodIndex){
		super(mass, 2, fragmentationMethodIndex);
		this.locationGroup = calculateMassGroup(mass, pm);
		this.setTol(tol);
	}
	
	@Override
	public boolean equals(Object obj) {
		if(obj == this) return true;
		if(obj instanceof Node){
			Node other = (Node)obj;
			return other.fragmentationMethodIndex == this.fragmentationMethodIndex && other.getMass() == this.getMass();
		}else return false;
	}
	
	@Override
	public float getAccuracy(){
		if(posteriorAccuracy >= 0) return posteriorAccuracy;
		return accuracy;
	}
	
	/**
	 * Gets the attachment cost.
	 *
	 * @return the attachment cost
	 */
	public float getAttachmentCost(){ // for prim's algorithm
		return attachmentCost;
	}
	
	/**
	 * Gets the FPV.
	 *
	 * @return the FPV
	 */
	public float[] getFPV() { return fpv; }


	/**
	 * Gets the location group.
	 *
	 * @return the location group
	 */
	public int getLocationGroup(){ return locationGroup; }
	
	@Override
	public float getLR(){
		if(posteriorLR > -Float.MAX_VALUE) return posteriorLR;
		return LR;
	}			
	
	/**
	 * Gets the min aa.
	 *
	 * @return the min aa
	 */
	public int getMinAA(){
		return minAA;
	}
	
	/**
	 * Gets the MS2 tol.
	 *
	 * @return the tol
	 */
	public Tolerance getTol() {
		return tol;
	}
	
	@Override
	public int hashCode() {
		return (int)(getMass()) * (fragmentationMethodIndex+1);
	}
	
	@Override
	public boolean isCorrect(Peptide pep, float offset) {
		boolean isCorrect = false;
		
		if(isSink){
			isCorrect = this.isWithinTolerance(pep.getMass() + offset, pep.getMass());
		}else{
			if(this.getMass() == offset) isCorrect =true;
			else{
				float prm = 0;
				for(AminoAcid aa : pep){
					prm += aa.getMass();
					if(this.isWithinTolerance(prm, pep.getMass())){
						isCorrect = true;
						break;
					}
				}
			}
		}
		return isCorrect;
	}
	
	@Override
	public boolean isEnzymatic() { return isEnzymaticNode; }
	
	/**
	 * Checks if is sink.
	 *
	 * @return true, if is sink
	 */
	public boolean isSink(){ return isSink; }
	
	/**
	 * Checks if is source.
	 *
	 * @return true, if is source
	 */
	public boolean isSource() { return this.getMass() == 0; }
	
	/**
	 * Sets the attachment cost.
	 *
	 * @param attachmentCost the new attachment cost
	 */
	public void setAttachmentCost(float attachmentCost){
		this.attachmentCost = attachmentCost;
	}
	
	@Override
	protected void setLRandAccuracy() {
		if(fpv !=null && ionWeights!=null && ionWeights[fragmentationMethodIndex] !=null && ionWeights[fragmentationMethodIndex][charge]!=null){
			int caseIndex=0;
			int cardinality = 0;
			int sigIonSize = fpv.length;			
			LR = 0;
			while(sigIonSize > 0){
				caseIndex=0;
				cardinality = 0;				
				for(int j=0; j<sigIonSize; j++){ // ion
					if(fpv[j] > 0){
						caseIndex += 1<<j;
						cardinality++;
						if(cardinality >= maxIonNumberPerNode) break;
					}
				}
				if(ionWeights[fragmentationMethodIndex][charge][caseIndex] != null) break;
				sigIonSize--;
			}			
			if(sigIonSize>0){
				int[] indices = new int[cardinality];
				int t = 0;				
				for(int j=0; j<sigIonSize; j++){ // ion
					if(fpv[j] > 0){
						indices[t++] = j;
						if(t >= indices.length) break;
					}
				}				
				accuracy = ionOffset[fragmentationMethodIndex][charge][caseIndex];				
				for(int j=0; j<cardinality; j++){
					float weighted = fpv[indices[j]] * ionWeights[fragmentationMethodIndex][charge][caseIndex][j];
					accuracy+= weighted;
				}
			}else if(isEnzymaticNode){
				isEmptyEnzymaticNode = true;
			}			
			nullAccuracy =  nodeNullAccuracies[fragmentationMethodIndex][charge][locationGroup];
			if(nodeNoPeakAccuracies!=null) noPeakAccuracy = nodeNoPeakAccuracies[fragmentationMethodIndex][charge][locationGroup];	
			float f = tol.getToleranceAsDa(500f)/tols[fragmentationMethodIndex].getToleranceAsDa(500f);
			nullAccuracy *= f;
			noPeakAccuracy *= f;
			accuracy = Math.min(1-nullAccuracy, accuracy);
			accuracy = Math.max(0, accuracy);
			if(this.getMass() == 0){
				accuracy = 1;
			}
			
			if(isEnzymaticNode){
				accuracy = Math.max(accuracy, Math.min(accuracy+1f/enzyme.getResidues().size()*0.9f, 0.9f));
			}
			
			LR = calculateLRfromAccuracy(nullAccuracy, accuracy);			
		}else if(fpv!=null){ // for training
			float j=1;
			for(float p : fpv) j *= 1-p;
			accuracy = LR = 1-j;
		}
	}
	
	/**
	 * Sets the min aa.
	 *
	 * @param a the new min aa
	 */
	public void setMinAA(int a){
		minAA = Math.min(minAA, a);
	}

	/**
	 * Sets the posterior lr and accuracy.
	 */
	public void setPosteriorLRnAccuracy(){
		if(this.isEmptyEnzymaticNode) return;
		
		if(getMinAA() < Integer.MAX_VALUE && nodePosteriorProbabilities != null && nodePosteriorProbabilities[fragmentationMethodIndex] != null && nodePosteriorProbabilities[fragmentationMethodIndex][charge] != null && nodePosteriorProbabilities[fragmentationMethodIndex][charge][0] != null){
			
			float p1 = nodePosteriorProbabilities[fragmentationMethodIndex][charge][0][Math.min(getMinAA(), nodePosteriorProbabilities[fragmentationMethodIndex][charge][0].length-1)];
			p1 *= getAccuracy();
			
			float p2 = nodePosteriorProbabilities[fragmentationMethodIndex][charge][1][Math.min(getMinAA(), nodePosteriorProbabilities[fragmentationMethodIndex][charge][1].length-1)];
			p2 *= 1-getAccuracy();
			
			p1 = Math.max(1e-6f, p1);
			p2 = Math.max(1e-6f, p2);
			
			posteriorAccuracy = p1/(p1+p2);
			posteriorLR = LR + (float)Math.log(p1/p2);
		}
	}



	/**
	 * Sets the sink.
	 *
	 * @param pmtol the new sink
	 */
	public void setSink(Tolerance pmtol){ isSink = true; tol = pmtol; this.accuracy = 1; LR = Math.min(100, calculateLRfromAccuracy(nullAccuracy, 1)); }

	/**
	 * Sets the MS2 tol.
	 *
	 * @param tol the new tol
	 */
	public void setTol(Tolerance tol) {
		this.tol = tol;
	}

	/**
	 * To file string.
	 *
	 * @return the string
	 */
	public String toFileString(){
		StringBuffer s = new StringBuffer();
		
		s.append(':');
		s.append(this.getMass());
		s.append(',');
		s.append(this.fragmentationMethodIndex);
		s.append(',');
		s.append(this.LR);
		s.append(':');
		
		return s.toString();
	}

	@Override
	public String toString(){
		StringBuffer s = new StringBuffer();
		
		s.append(':');
		s.append(this.getMass());
		s.append(',');
		s.append(this.fragmentationMethodIndex);
		s.append(',');
		s.append(this.LR);
		s.append(':');
		
		return s.toString();
	}
}
