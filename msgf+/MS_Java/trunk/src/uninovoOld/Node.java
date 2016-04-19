package uninovoOld;



import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.Peptide;



public class Node extends SpectrumGraphComponent{
	public final static int maxCardinality = 3;
	
	private boolean isEnzymaticNode = false;
	protected boolean isEmptyEnzymaticNode = false;
	private boolean isSink = false;
	private float[] ionProbs;
	//private float[] noiseProbs;
	private float attachmentCost = Float.MAX_VALUE;
	private int caseNumber = 0;
	private int section = 0;
	private int minAA = Integer.MAX_VALUE;
	private float posteriorLR=-Float.MAX_VALUE;
	private float posteriorAccuracy = -1;
	//private float[] ionProbsSelected;
	
	public Node(float mass, int charge, float[] ionProbs, Tolerance tol, float pm, boolean isEnzymaticNode, int typeIndex){
		super(mass, charge, typeIndex);
		this.ionProbs = ionProbs;		
		this.isEnzymaticNode = isEnzymaticNode;
		this.section = calculateSection(mass, pm);
		this.setTol(tol);
		setLRandAccuracy();
	}
	
	public Node(float mass, Tolerance tol, float pm, int typeIndex){
		super(mass, 2, typeIndex);
		this.section = calculateSection(mass, pm);
		this.setTol(tol);
	}
	
	public int getSection(){ return section; }
	
	/*protected void updatePosteriorAccuracy(float a){ 
		if(this.posteriorAccuracy < a){
			this.posteriorAccuracy = a;
			//this.posteriorLR = calulateLRfromAccuracy(nullAccuracy, posteriorAccuracy, false);
		}
	}*/
	

	
	public int getMinAA(){
		return minAA;
	}
	
	public void setMinAA(int a){
		minAA = Math.min(minAA, a);
	}
	
	@Override
	public float getAccuracy(){
		if(posteriorAccuracy >= 0) return posteriorAccuracy;
		return accuracy;
	}
	
	@Override
	public float getLR(){
		if(posteriorLR > -Float.MAX_VALUE) return posteriorLR;
		return LR;
	}
	
	public void setPosteriorLRnAccuracy(){
		if(this.isEmptyEnzymaticNode) return;
		
		if(getMinAA() < Integer.MAX_VALUE && nodePosteriorProbabilities != null && nodePosteriorProbabilities[typeIndex] != null && nodePosteriorProbabilities[typeIndex][charge] != null && nodePosteriorProbabilities[typeIndex][charge][0] != null){
			
			float p1 = nodePosteriorProbabilities[typeIndex][charge][0][Math.min(getMinAA(), nodePosteriorProbabilities[typeIndex][charge][0].length-1)];
			p1 *= getAccuracy();
			
			float p2 = nodePosteriorProbabilities[typeIndex][charge][1][Math.min(getMinAA(), nodePosteriorProbabilities[typeIndex][charge][1].length-1)];
			p2 *= 1-getAccuracy();
			
			p1 = Math.max(1e-6f, p1);
			p2 = Math.max(1e-6f, p2);
			
			posteriorAccuracy = p1/(p1+p2);
			//System.out.println(posteriorAccuracy);
		//	System.out.println(calulateLRfromAccuracy(nullAccuracy, posteriorAccuracy, false));
			posteriorLR = LR + (float)Math.log(p1/p2);
		//	posteriorLR = calulateLRfromAccuracy(nullAccuracy, posteriorAccuracy, false);
		//	System.out.println(posteriorLR+"*");
			//updatePosteriorLR(Math.max(1e-6f, p1) / (p1+p2));
		}
	}

	/*public Node(float mass, int charge, float LR, int typeIndex){ // just for the ROC test. erase later...
		super(mass, charge, typeIndex);
		this.LR = LR;
	}*/
	/*
	private void selectSigIonProbs(int numSig){
		int n=0;
		int[] numPrefix = new int[this.charge];
		int[] numSuffix = new int[this.charge];
		float[] ionProbsSelected = new float[numSig];
		int caseNumber = 0;
		
		if(ionProbs != null){
			
			for(int k=0; k< ionProbs.length; k++){
				float i = ionProbs[k];
				
				if(i>0){
					int c = ions.get(k).getCharge();
					int[] num;
					if(ions.get(k) instanceof IonType.PrefixIon){
						num = numPrefix;
					}else{
						num = numSuffix;
					}
					
					isNoPeakNode = false;
					
					if(num[c] >= 1){
						continue;
					}
					
					if(n>=numSig){
						continue;
					}
					
					
					ionProbsSelected[n] = i;
					caseNumber += 1<<k;
					
					
					if(ionProbProbabilities!= null && ionProbProbabilities[charge]!=null && n>0){
						int length = ionProbProbabilities[charge][caseNumber].length;
						int bin1 = Math.round((length-1)*ionProbsSelected[0]); 
						int bin2 = Math.round((length-1)*ionProbsSelected[1]); 
						
						if(!isEnzymaticNode){
							if(ionProbProbabilities[charge][caseNumber][bin1][bin2][1] == 0){
								ionProbsSelected[n] = 0;
								caseNumber -= 1<<k;
							}
						}else{
							if(enzymaticIonProbProbabilities[charge][caseNumber][bin1][bin2][1] == 0){
								ionProbsSelected[n] = 0;
								caseNumber -= 1<<k;
							}
						}
					}

					num[c]++; n++;
				}
			}
			this.ionProbsSelected = ionProbsSelected;
			this.caseNumber = caseNumber;
		}
	}
	*/
	public float[] getIonProbs() { return ionProbs; }
	
//	public float[] getIonProbsSelected() { return ionProbsSelected; }
	
	public int getCaseNumber() { return caseNumber; }

	public void setSink(Tolerance pmtol){ isSink = true; tol = pmtol; this.accuracy = 1; LR = Math.min(100, calulateLRfromAccuracy(nullAccuracy, 1, false)); }			
	
	public boolean isSink(){ return isSink; }
	
	public boolean isSource() { return this.getMass() == 0; }
	
	@Override
	public boolean isEnzymatic() { return isEnzymaticNode; }
	
	@Override
	public String toString(){
		StringBuffer s = new StringBuffer();
		
		s.append(':');
		s.append(this.getMass());
		s.append(',');
		s.append(this.typeIndex);
		s.append(',');
		s.append(this.LR);
		s.append(':');
		
		//s += this.getMass()+" ";//+"\t"+this.getAccuracy();
		return s.toString();
	}
	
	public String toFileString(){
		StringBuffer s = new StringBuffer();
		
		s.append(':');
		s.append(this.getMass());
		s.append(',');
		s.append(this.typeIndex);
		s.append(',');
		s.append(this.LR);
		s.append(':');
		
		//s += this.getMass()+" ";//+"\t"+this.getAccuracy();
		return s.toString();
	}
	
	@Override
	protected void setLRandAccuracy() {
		
		if(ionProbs !=null && ionWeights!=null && ionWeights[typeIndex] !=null && ionWeights[typeIndex][charge]!=null){
			int caseIndex=0;
			int cardinality = 0;
			int sigIonSize = ionProbs.length;			
			//float maxIonProb = 0;
			LR = 0;
			while(sigIonSize > 0){
				caseIndex=0;
				cardinality = 0;
				
				for(int j=0; j<sigIonSize; j++){ // ion
					if(ionProbs[j] > 0){
						caseIndex += 1<<j;
						cardinality++;
						if(cardinality >= maxCardinality) break;
					}
				}
				if(ionWeights[typeIndex][charge][caseIndex] != null) break;
				sigIonSize--;
			}
			
			if(sigIonSize>0){
				int[] indices = new int[cardinality];
				int t = 0;
				
				for(int j=0; j<sigIonSize; j++){ // ion
					if(ionProbs[j] > 0){
					//	maxIonProb = Math.max(maxIonProb, ionProbs[j]);
						indices[t++] = j;
						if(t >= indices.length) break;
					}
				}
				
				accuracy = ionOffset[typeIndex][charge][caseIndex];
				
				for(int j=0; j<cardinality; j++){
					float weighted = ionProbs[indices[j]] * ionWeights[typeIndex][charge][caseIndex][j];
					accuracy+= weighted;
				//	LR += Math.log(Math.max(1e-3f, weighted)/Math.max(1e-3f, noiseProbs[indices[j]]));
				}
				//accuracy = Math.max(accuracy, 0.01f);
				
			}else if(isEnzymaticNode){
				isEmptyEnzymaticNode = true;
			}
			
			nullAccuracy =  nodeNullAccuracies[typeIndex][charge][section];//2f/121.6f*tol.getToleranceAsDa(500f);
			if(nodeNoPeakAccuracies!=null) noPeakAccuracy = nodeNoPeakAccuracies[typeIndex][charge][section];
			
			float f = tol.getToleranceAsDa(500f)/tols[typeIndex].getToleranceAsDa(500f);
			nullAccuracy *= f;
			noPeakAccuracy *= f;
		//	System.out.println(sectionNum + "\t" + nullAccuracy + "\t" + noPeakAccuracy);
			accuracy = Math.min(1-nullAccuracy, accuracy);
			accuracy = Math.max(0, accuracy);
			
			if(this.getMass() == 0){
				accuracy = 1;
			}
			
			//nullAccuracy = nodesNullAccuracies[typeIndex][charge];// 1f/121f;//
				
			if(isEnzymaticNode){
				accuracy = Math.max(accuracy, Math.min(accuracy+1f/enzyme.getResidues().size()*0.9f, 0.9f));//TODO .. do something
			}
			
			LR = calulateLRfromAccuracy(nullAccuracy, accuracy, false);
			//System.out.println(LR);
			//noiseProbs
		
	/*		for(int i=0;i<ionProbs.length;i++){
				if(isEnzymaticNode){
					LR += Math.log(0.8f/Math.max(1e-1f, noiseProbs[i]));
				}else{
					LR += Math.log(Math.max(1e-1f, ionProbs[i])/Math.max(1e-1f, noiseProbs[i]));
					//System.out.println( ionProbs[i] + "\t" +  noiseProbs[i]);
				}
			}
		*/
		//	if(!this.isSink && this.getMass() !=0 && accuracy > 0.7f) System.out.println("* " + this + "\t" + 
				//	"\t" + accuracy + "\t" + LR);
			
			//accuracy = Math.min(accuracy, 0.97f);
			
			
			
		}else if(ionProbs!=null){ // for training
			float j=1;
			for(float p : ionProbs) j *= 1-p;
			accuracy = LR = 1-j;
		}
	}
	
	static private float calulateLRfromAccuracy(double nacc, double acc, boolean isNoPeak){
		double up = 1-nacc;//, down = 1e-1;
		
		if(isNoPeak)
			acc = getBoundedValue(acc, up, 0);
		else
			acc = getBoundedValue(acc, up, nacc);
		
		return (float) (Math.log(acc/(1-acc) * (1-nacc)/nacc));
		
	}

	static public float getNoPeakLR(int typeIndex, int charge, float mass, float pm){
		charge = Math.min(Math.max(chargeRanges[typeIndex][0] , charge), chargeRanges[typeIndex][1]);
		int section = calculateSection(mass, pm);
		if(nodeNoPeakAccuracies == null || nodeNoPeakAccuracies[typeIndex]==null || nodeNoPeakAccuracies[typeIndex][charge] == null || nodeNoPeakAccuracies[typeIndex][charge][section] == 0){
			return calulateLRfromAccuracy(0.1, 0.1/3, true);
		}else{
			float n = nodeNoPeakAccuracies[typeIndex][charge][section];
			return calulateLRfromAccuracy(nodeNullAccuracies[typeIndex][charge][section],  n, true);
		}

	}
	
	static private int calculateSection(float mass, float pm){
		return Math.max(0, Math.min(sectionNumber-1, (int)(mass/pm * sectionNumber)));
	}
	
	public float getAttachmentCost(){ // for prim's algorithm
		return attachmentCost;
	}
	
	public void setAttachmentCost(float attachmentCost){
		this.attachmentCost = attachmentCost;
	}

	@Override
	public boolean equals(Object obj) {
		if(obj == this) return true;
		if(obj instanceof Node){
			Node other = (Node)obj;
			return other.typeIndex == this.typeIndex && other.getMass() == this.getMass();
		}else return false;
	}



	@Override
	public int hashCode() {
		return (int)(getMass()) * (typeIndex+1);//Float(getMass()).hashCode();// * new Integer(typeIndex).hashCode();
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

	public void setTol(Tolerance tol) {
		this.tol = tol;
	}

	public Tolerance getTol() {
		return tol;
	}
}
