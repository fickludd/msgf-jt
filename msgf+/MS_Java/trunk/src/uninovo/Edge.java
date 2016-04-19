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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import uninovo.util.AminoAcid;
import uninovo.util.AminoAcidSet;
import uninovo.util.Constants;
import uninovo.util.ModifiedAminoAcid;
import uninovo.util.NominalMass;
import uninovo.util.Peptide;

/**
 * The Class Edge represents an edge in a de novo sequence
 */
public class Edge  extends SpectrumGraphComponent{
	
	/** The aa composition table. */
	private static HashMap<Double, HashSet<Character>> aaCompTable;
	
	/** The aa residue table. */
	private static HashMap<Double, AminoAcid> aaResidueTable;
	
	/** The amino acid set. */
	static private AminoAcidSet aaSet;
	
	/** The Constant massDeviationIndexNum. */
	static final public int massDeviationIndexNum = 80;
	
	/** The max aa num table. */
	private static HashMap<Double, Integer> maxAANumTable;
	
	/** The Constant max edge mass, i.e., max gap mass. */
	static private final float maxEdgeMass = 900;//600
	
	/** The min aa num table. */
	private static HashMap<Double, Integer> minAANumTable;
	
	/** The min edge mass. */
	static private float minMass = 57;
	
	/** The modification table. */
	private static HashMap<Integer, ArrayList<ArrayList<AminoAcid>>> modTable;
	
	/** The possible molecular masses (consisting of amino acids). */
	private static ArrayList<Double> possibleMolecularMasses; // sorted in ascending order
	
	/**
	 * Update necessary tables.
	 *
	 * @param aaSet the amino acid set
	 */
	private static void updateTables(AminoAcidSet aaSet){ // revise.. sometime!
		Edge.aaSet = aaSet;
	
		possibleMolecularMasses = new ArrayList<Double>();
		minAANumTable = new HashMap<Double, Integer>();
		maxAANumTable = new HashMap<Double, Integer>();
		aaResidueTable = new HashMap<Double, AminoAcid>();
		aaCompTable = new HashMap<Double, HashSet<Character>>();
		modTable = new HashMap<Integer, ArrayList<ArrayList<AminoAcid>>>();
		HashMap<Integer, ArrayList<AminoAcid>> modAAs = new HashMap<Integer, ArrayList<AminoAcid>>();
		ArrayList<Integer> possibleUnModSeqMasses = new ArrayList<Integer>();
		HashSet<Integer> unModAAMasses = new HashSet<Integer>();
		
		for(AminoAcid aa : aaSet){
			double m = aa.getAccurateMass();
			m =  (Math.round(m*Constants.INTEGER_MASS_SCALER_HIGH_PRECISION))/Constants.INTEGER_MASS_SCALER_HIGH_PRECISION;

			aaResidueTable.put(m, aa);
			
			if(aa.isModified()){
				ArrayList<AminoAcid> aas = new ArrayList<AminoAcid>();
				int nm = aa.getNominalMass();
				if(!modAAs.containsKey(nm))
					modAAs.put(nm,  new ArrayList<AminoAcid>());
				aas = modAAs.get(nm);
				aas.add(aa);
				
			}else{
				unModAAMasses.add(aa.getNominalMass());
			}
			
		}

		possibleUnModSeqMasses.add(0);
		
		for(int i=0; i< possibleUnModSeqMasses.size(); i++){
			int m = possibleUnModSeqMasses.get(i);
			for(int aaMass : unModAAMasses){
				int mass = m + aaMass;
				int k = Collections.binarySearch(possibleUnModSeqMasses, mass);
				if(k<0){
					possibleUnModSeqMasses.add( -k-1,mass);
				}
			}
			
			Collections.sort(possibleUnModSeqMasses);
			if(possibleUnModSeqMasses.get(possibleUnModSeqMasses.size()-1) > maxEdgeMass) break;
		}
			
		minAANumTable.put(0.0, 0);
		maxAANumTable.put(0.0, 0);
		possibleMolecularMasses.add(0.0);
		
		for(int i=0; i< possibleMolecularMasses.size(); i++){
			double m = possibleMolecularMasses.get(i);
			for(double aaMass : aaResidueTable.keySet()){
				double mass = m + aaMass;
			
				int k = Collections.binarySearch(possibleMolecularMasses, mass);
				if(k<0){
					possibleMolecularMasses.add( -k-1,mass);
				}
		
				HashSet<Character> aaPrevComp = null;
				
				if(aaCompTable.containsKey(m))
					aaPrevComp = new HashSet<Character>(aaCompTable.get(m));
				else 
					aaPrevComp = new HashSet<Character>();
				
				aaPrevComp.add(aaResidueTable.get(aaMass).getResidue());		
				if(aaCompTable.containsKey(mass)) {
					aaPrevComp.addAll(aaCompTable.get(mass));
				}
				aaCompTable.put(mass, aaPrevComp);
				
				if(modAAs.containsKey(NominalMass.toNominalMass((float)aaMass))){
					ArrayList<ArrayList<AminoAcid>> prevMod = null;
					int nm = NominalMass.toNominalMass((float)m);
					int nmass = NominalMass.toNominalMass((float)mass);
					
					if(modTable.containsKey(nm))
						prevMod = modTable.get(nm);
					else{
						prevMod = new ArrayList<ArrayList<AminoAcid>>();
						prevMod.add(new ArrayList<AminoAcid>());
					}
					
					if(!modTable.containsKey(nmass)){
						modTable.put(nmass, new ArrayList<ArrayList<AminoAcid>>());
						if(Collections.binarySearch(possibleUnModSeqMasses, nmass)>=0) 
							modTable.get(nmass).add(new ArrayList<AminoAcid>());
					}
					
					ArrayList<ArrayList<AminoAcid>> currMod = modTable.get(nmass);
					
					for(AminoAcid maa : modAAs.get(NominalMass.toNominalMass((float)aaMass))){
						for(ArrayList<AminoAcid> aas : prevMod){
							ArrayList<AminoAcid> naas = new ArrayList<AminoAcid>(aas);
							naas.add(maa);
							Collections.sort(naas);
							
							if(!currMod.contains(naas) && naas.size() <= aaSet.getMaxNumberOfVariableModificationsPerPeptide())
								currMod.add(naas);
						}
					}
					
				}	
				
				int minAANum = minAANumTable.get(m) + 1;
				int maxAANum = maxAANumTable.get(m) + 1;
				if(minAANumTable.containsKey(mass)){
					minAANum = Math.min(minAANum, minAANumTable.get(mass));
				}
				if(maxAANumTable.containsKey(mass)){
					maxAANum = Math.max(maxAANum, maxAANumTable.get(mass));
				}
				minAANumTable.put(mass, minAANum);
				maxAANumTable.put(mass, maxAANum);
			
			}
			
			if(possibleMolecularMasses.get(possibleMolecularMasses.size()-1) > maxEdgeMass) break;
		}
			
		Collections.sort(possibleMolecularMasses);
	}
	
	/** The closest molecular mass. */
	private Double closestMolecularMass=0.0;
	
	/** Is this exceed the max mass? */
	private boolean exceedsMaxMass = false;
	
	/** Is this a composite mass? */
	private boolean isCompositeMass = false;
	
	/** Is this enzymatic edge? */
	private boolean isEnzymaticEdge = false;
	
	/** Is this edge gap? */
	private boolean isGap = true;
	
	/** Is this a PTM edge?. */
	private boolean isPTM = false;
	
	/** Is this a valid edge? */
	private boolean isValid = true;
	
	/** The left and right nodes */
	private Node l, r;
	
	/** The mass deviation index. */
	private int massDeviationIndex = Integer.MIN_VALUE;
	
	/** The min and max amino acid number that can have the mass of this edge */
	private int minAANum, maxAANum;

	/**
	 * Instantiates a new dummy edge. Only used in the dynamic programming to generate sequences
	 *
	 * @param typeIndex the type index
	 */
	public Edge(int typeIndex){
		super(0, 1, typeIndex);
		isValid = false;
	}
	
	
	/**
	 * Instantiates a new edge.
	 *
	 * @param l the left node
	 * @param r the right node
	 * @param pm the parent mass
	 * @param aaSet the amino acid set
	 */
	public Edge(Node l, Node r, float pm, AminoAcidSet aaSet) {
		super(r.getMass() - l.getMass(), l.charge, l.fragmentationMethodIndex);
	
		this.l = l;
		this.r = r;
		
		this.tol = l.getTol();
		
		if(this.tol.getToleranceAsDa(100) < r.getTol().getToleranceAsDa(100)){
			this.tol = r.getTol();
			this.fragmentationMethodIndex = r.fragmentationMethodIndex;
			this.charge = r.charge;
		}

		
		if(Edge.aaSet != aaSet) updateTables(aaSet);
		
		if(getMass() > maxEdgeMass){
			exceedsMaxMass = true;
		}
		
		// gather masses within Tolerance tol
		float t = tol.getToleranceAsDa(Math.max(pm - l.getMass(), r.getMass()));
		t = Math.min(t, 0.5f);
		
		if(!exceedsMaxMass){
			int i1 = Collections.binarySearch(possibleMolecularMasses, (double)(getMass()-t));
			int i2 = Collections.binarySearch(possibleMolecularMasses, (double)(getMass()+t));
		
			i1 = i1<0? -i1-1 : i1;
			i2 = i2<0? -i2-1 : i2;
			
			minAANum = Integer.MAX_VALUE;
			maxAANum = Integer.MIN_VALUE;
			
			HashSet<Character> comp = new HashSet<Character>();
			
			for(int i=Math.max(0, i1-1);i<Math.min(i2+1, possibleMolecularMasses.size());i++){
				double m = possibleMolecularMasses.get(i);
				
				float d = Math.abs((float)m - getMass());
				if(d<1.5f*t){
					int minAANumt = minAANumTable.get(m);
					HashSet<Character> subcomp = aaCompTable.get(m);
					if(subcomp != null) comp.addAll(subcomp);
					
					if(minAANumt < minAANum){
						minAANum = Math.min(minAANum, minAANumt);
						closestMolecularMass = m;
					}else if(minAANumt == minAANum && d<Math.abs(closestMolecularMass - getMass())){
						closestMolecularMass = m;
					}
					maxAANum = Math.max(maxAANum, maxAANumTable.get(m));
				}
				
			}
			isGap = minAANum > 1;
			isCompositeMass = (!isGap) && (maxAANum > 1);
			
			if(enzyme != null){
				if(exceedsMaxMass) isEnzymaticEdge = true;
				else if((enzyme.isCTerm() && r.isSink()) || (enzyme.isNTerm() && l.getMass() == 0)){
					for(AminoAcid aa : enzyme.getResidues()){
						if(comp.contains(aa.getResidue())){
							isEnzymaticEdge = true;
							break;
						}
					}
				}
			}
		}else{
			minAANum = Integer.MAX_VALUE;
			maxAANum = Integer.MAX_VALUE;
			closestMolecularMass = (double)getMass();
			isGap = true;
			isCompositeMass = false;
		}
		
		float diff = (float)(closestMolecularMass - getMass());

		if(closestMolecularMass == 0 || Math.abs(diff) > 1.5f*t){
			//isValid = false;
			if(!exceedsMaxMass) isPTM = true;
			minAANum = Integer.MAX_VALUE;
			maxAANum = Integer.MAX_VALUE;
			closestMolecularMass = (double)getMass();
			isGap = true;
			isCompositeMass = false;
		}
	
		if(getMass() < minMass) {
			isValid = false;
		}
		
		if(minAANum< Integer.MAX_VALUE && !SpectrumGraphComponent.useBinning(tol)){
			massDeviationIndex = Math.round(diff/(1.5f*t)*massDeviationIndexNum/2+massDeviationIndexNum/2);
			massDeviationIndex = Math.max(0, massDeviationIndex);
			massDeviationIndex = Math.min(massDeviationIndexNum-1, massDeviationIndex);
		}
	}

	@Override
	public boolean equals(Object obj) {
		if(obj == this) return true;
		if(obj instanceof Edge){
			Edge o = (Edge)obj;			
			return o.l.equals(this.l) && o.r.equals(this.r);
		}else return false;
	}
	
	/**
	 * Gets the left node.
	 *
	 * @return the left node
	 */
	public Node getLeftNode() { return l; }
	
	@Override
	public float getLR(){
		return this.r.getLR() + LR;
	}
	
	@Override
	public float getMass(){
		float  mass = super.getMass();
		if(this.l.isSource()){
			for(ModifiedAminoAcid aa : aaSet.getNTermFixedMods()){
				mass -= aa.getModification().getMass();
			}
		}
		if(this.r.isSink()){
			for(ModifiedAminoAcid aa : aaSet.getCTermFixedMods()){
				mass -= aa.getModification().getMass();
			}
		}
		return Math.max(0, mass);
	}
	
	/**
	 * Gets the mass deviation index - parameters are separately learned for edges with different mass deviation indices
	 *
	 * @return the mass deviation index
	 */
	public int getMassDeviationIndex() { return massDeviationIndex; }
	
	/**
	 * Gets the min aa num.
	 *
	 * @return the min aa num
	 */
	public int getMinAANum() { return minAANum; }
	
	/**
	 * Gets the modified amino acid list, if there is.
	 *
	 * @return the modified aa list
	 */
	public ArrayList<ArrayList<AminoAcid>> getModifiedAAList(){
		return modTable.get(NominalMass.toNominalMass(this.getMass()));
	}
	
	/**
	 * Gets the right node.
	 *
	 * @return the right node
	 */
	public Node getRightNode() { return r; }
	
	@Override
	public int hashCode() {
		if(l==null || r==null) return 0;
		return l.hashCode() + r.hashCode();
	}
	
	/**
	 * Checks if this is composite mass.
	 *
	 * @return true, if this is composite mass
	 */
	public boolean isCompositeMass() { return isCompositeMass; }
	
	@Override
	public boolean isCorrect(Peptide pep, float offset){
		boolean isCorrect = this.l.isCorrect(pep, offset);
		isCorrect &= this.r.isCorrect(pep, offset);
		
		if(isCorrect && isCompositeMass){
			float prm = 0;
			for(AminoAcid aa : pep){
				prm += aa.getMass();
				if(prm < l.getMass() + 1) continue;
				if(prm > r.getMass() - 1) break;
				
				isCorrect = false;
				break;
				
			}
		}
	
		return isCorrect;
		
	}
	
	@Override
	public boolean isEnzymatic() { return isEnzymaticEdge; }

	/**
	 * Checks if this is gap.
	 *
	 * @return true, if this is gap
	 */
	public boolean isGap() { return isGap; }

	/**
	 * Checks if this is PTM.
	 *
	 * @return true, if this is PTM
	 */
	public boolean isPTM() { return isPTM; }
	
	/**
	 * Checks if this is valid.
	 *
	 * @return true, if this is valid
	 */
	public boolean isValid() { return isValid; }
	
	@Override
	public void setLRandAccuracy() {
		float acc;
		float accR = r.getAccuracy();
		float accL = l.getAccuracy();
	
		int i1 = minAANum;
		
		if(r.isSink()){
			acc = accL; 
		}else if(l.getMass() == 0){
			acc = accR;
		}else if(edgeAccuracies != null && edgeAccuracies[fragmentationMethodIndex] != null && edgeAccuracies[fragmentationMethodIndex][charge] != null){
			 if(!isValid) acc = 0;
			 else{
				i1 = Math.min(edgeAccuracies[fragmentationMethodIndex][charge].length-1, i1);
				if(isCompositeMass) i1 = 0;
				int i2_1 = (int) ((edgeAccuracies[fragmentationMethodIndex][charge][0].length-1)*accR);
				int i2_0 = i2_1-1;
				int i2_2 = i2_1 + 1;
				int i3_1 = (int)((edgeAccuracies[fragmentationMethodIndex][charge][0].length-1)*accL);
				int i3_2 = i3_1 + 1;
				if(i3_2 > edgeAccuracies[fragmentationMethodIndex][charge][0].length-1) i3_2 = i3_1-1;
				
				
				float acct0 =0 , acct1 = 0; 

				float acc1 = edgeAccuracies[fragmentationMethodIndex][charge][i1][i2_1][i3_1][0];
				float acc2 = 0, acc0 = 0;
				
				if(i2_0<0){
					acc0 = 0;
				}else
					acc0 = edgeAccuracies[fragmentationMethodIndex][charge][i1][i2_0][i3_1][0];
				
				if(i2_2 >= edgeAccuracies[fragmentationMethodIndex][charge][i1].length){
					acc2 = accR;
				}else
					acc2 = edgeAccuracies[fragmentationMethodIndex][charge][i1][i2_2][i3_1][0];
				
				acct0 = acc1 + (acc2-acc1)/(i2_2-i2_1)*((edgeAccuracies[fragmentationMethodIndex][charge][0].length-1)*accR - i2_1);
				acct0 += acc1 + (acc0-acc1)/(i2_0-i2_1)*((edgeAccuracies[fragmentationMethodIndex][charge][0].length-1)*accR - i2_1);
				acct0/=2;
				
				acc1 = edgeAccuracies[fragmentationMethodIndex][charge][i1][i2_1][i3_2][0];
				acc2 = 0; acc0 = 0;
	
				if(i2_0<0){
					acc0 = 0;
				}else
					acc0 = edgeAccuracies[fragmentationMethodIndex][charge][i1][i2_0][i3_2][0];
				
				if(i2_2 >= edgeAccuracies[fragmentationMethodIndex][charge][i1].length){
					acc2 = accR;
				}else
					acc2 = edgeAccuracies[fragmentationMethodIndex][charge][i1][i2_2][i3_2][0];
			
				acct1 = acc1 + (acc2-acc1)/(i2_2-i2_1)*((edgeAccuracies[fragmentationMethodIndex][charge][0].length-1)*accR - i2_1);
				acct1 += acc1 + (acc0-acc1)/(i2_0-i2_1)*((edgeAccuracies[fragmentationMethodIndex][charge][0].length-1)*accR - i2_1);
				acct1/=2;
				
				acc = acct0 + (acct1-acct0)/(i3_2-i3_1)*((edgeAccuracies[fragmentationMethodIndex][charge][0].length-1)*accL - i3_1);
				

			 }
		
		}else{
			acc = accR * accL;// nacc = 0;
		}
		
		LR = 0;
		
		acc = Math.min(acc, 1-r.nullAccuracy);// for stability
		acc = Math.max(acc, 0);

		
		if(!this.r.isEmptyEnzymaticNode && !this.l.isEmptyEnzymaticNode){
			if(edgeDeviationAccuracies !=null && edgeDeviationAccuracies[fragmentationMethodIndex]!=null && edgeDeviationAccuracies[fragmentationMethodIndex][charge]!=null && edgeDeviationAccuracies[fragmentationMethodIndex][charge][0]!=null  && massDeviationIndex >=0){
				float pp1 = edgeDeviationAccuracies[fragmentationMethodIndex][charge][0][Math.min(edgeDeviationAccuracies[fragmentationMethodIndex][charge][0].length-1, minAANum)][massDeviationIndex];
				float pp2 = edgeDeviationAccuracies[fragmentationMethodIndex][charge][1][Math.min(edgeDeviationAccuracies[fragmentationMethodIndex][charge][0].length-1, minAANum)][massDeviationIndex];
				
				acc = pp1*acc/(pp1*acc + pp2*(1-acc));
				LR += Math.log(pp1/pp2);
			}
		}
		
		if(isPTM()) LR += Math.log(0.05f);
		
		accuracy = Math.max(acc, accL*accR);
		
		float factor = 0.3f;
		
		accuracy = Math.min(accuracy, accR + (float)Math.sqrt(accR * (1-accR)) * factor);//
		accuracy = Math.min(accuracy, accL + (float)Math.sqrt(accL * (1-accL)) * factor);// 
		
	}
	
	@Override
	public String toString(){
		StringBuffer s = new StringBuffer();
		if(isGap){
			s.append('[');
			s.append(getMass());
			if(isPTM) s.append('*');
			s.append(']');
		}else{
			s.append(aaResidueTable.get(closestMolecularMass).getResidueStr());
		}
		return s.toString();
	}
	
	/**
	 * Weight for minimum spanning tree. Used to calculate accuracy of sequences
	 *
	 * @return the float
	 */
	public float weightForMST(){
		float s = this.getAccuracy();//
		float t = l.getAccuracy() + r.getAccuracy() - s;
		t = Math.max(t, 0);
		t = Math.min(t, 1);
		return t;
	}

}
