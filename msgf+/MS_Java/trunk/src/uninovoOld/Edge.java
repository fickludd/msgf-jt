package uninovoOld;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import msgf.NominalMass;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Constants;
import msutil.ModifiedAminoAcid;
import msutil.Peptide;

//TODO possibleMolMasses, possibleNMolMasses, possibleCMolMasses, aaResidueTable, aaResidueNTable, C, modTable, modNTable, C

public class Edge  extends SpectrumGraphComponent{
	static private AminoAcidSet aaSet;
	static private final float maxEdgeMass = 900;//600
	static private float minMass = 57;

	//static public final float maxPepMass = 1200;
	static final public int massDeviationIndexNum = 80;
	
	private boolean isGap = true;
	private boolean isPTM = false;
	private boolean isCompositeMass = false;
	private boolean isValid = true;
	private boolean exceedsMaxMass = false;
	private boolean isEnzymaticEdge = false;
	private Node l, r;
	private Double closestMolecularMass=0.0;
	private int massDeviationIndex = Integer.MIN_VALUE;
	private int minAANum, maxAANum;
	
	private static ArrayList<Double> possibleMolecularMasses; // sorted in ascending order
	private static HashMap<Double, Integer> minAANumTable;
	private static HashMap<Double, Integer> maxAANumTable;
	private static HashMap<Double, AminoAcid> aaResidueTable;
	private static HashMap<Double, HashSet<Character>> aaCompTable;
	private static HashMap<Integer, ArrayList<ArrayList<AminoAcid>>> modTable;
	
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
				//mass =  (Math.round(mass*Constants.INTEGER_MASS_SCALER_HIGH_PRECISION))/Constants.INTEGER_MASS_SCALER_HIGH_PRECISION;

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

	public Edge(Node l, Node r, float pm, AminoAcidSet aaSet) {
		super(r.getMass() - l.getMass(), l.charge, l.typeIndex);
	
		this.l = l;
		this.r = r;
		
		this.tol = l.getTol();
		
		if(this.tol.getToleranceAsDa(100) < r.getTol().getToleranceAsDa(100)){
			this.tol = r.getTol();
			this.typeIndex = r.typeIndex;//0.9366197	133	142	8.819549
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
			
			minAANum = Integer.MAX_VALUE;//minAANumTable.get(closestMolecularMass);
			maxAANum = Integer.MIN_VALUE;//maxAANumTable.get(closestMolecularMass);
			
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
	//	if(isGap && (getLeftNode().getMass() !=0 &&!getRightNode().isSink()) ){
	//		//isValid = false;//tmp TODO
	//	}//else System.out.println(this);
		
		//massDeviationIndex = (int) Math.min(massDeviationIndexNum-1, (int)(diff/(t*1.5f) * massDeviationIndexNum));
		
	//	if(!isGap && this.getMass()>0) System.out.println(this + "\t" + closestMolecularMass + "\t" + minAANum + "\t" + maxAASeqNum);
		//setLRandAccuracy();
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
	
	
	public Edge(int typeIndex){
		super(0, 1, typeIndex);
		isValid = false;
	}
	
	@Override
	public boolean isEnzymatic() { return isEnzymaticEdge; }
	
	public Node getLeftNode() { return l; }
	public Node getRightNode() { return r; }
	public boolean isValid() { return isValid; }
	//public boolean exceedsMaxMass() {return exceedsMaxMass; }
	public boolean isGap() { return isGap; }
	public boolean isPTM() { return isPTM; }
	public boolean isCompositeMass() { return isCompositeMass; }
	public int getMassDeviationIndex() { return massDeviationIndex; }
	public float weightForMST(){
		float s = this.getAccuracy();//
		float t = l.getAccuracy() + r.getAccuracy() - s;
		t = Math.max(t, 0);
		t = Math.min(t, 1);
		return t;
	}
	public int getMinAANum() { return minAANum; }
	
	@Override
	public boolean equals(Object obj) {
		if(obj == this) return true;
		if(obj instanceof Edge){
			Edge o = (Edge)obj;			
			return o.l.equals(this.l) && o.r.equals(this.r);
		}else return false;
	}

	@Override
	public int hashCode() {
		if(l==null || r==null) return 0;
		return l.hashCode() + r.hashCode();
	}

	@Override
	public void setLRandAccuracy() {
		float acc;//, nacc = 0;
		float accR = r.getAccuracy();
		float accL = l.getAccuracy();
	
		int i1 = minAANum;
		
		if(r.isSink()){
			acc = accL; //nacc = 0;
		}else if(l.getMass() == 0){
			acc = accR;// nacc = 0;
		//}else if(exceedsMaxMass){
			//acc = accR * accL;			
		}else if(edgeAccuracies != null && edgeAccuracies[typeIndex] != null && edgeAccuracies[typeIndex][charge] != null){
			 if(!isValid) acc = 0;
			 else{
				i1 = Math.min(edgeAccuracies[typeIndex][charge].length-1, i1);
				if(isCompositeMass) i1 = 0;
				int i2_1 = (int) ((edgeAccuracies[typeIndex][charge][0].length-1)*accR);
				int i2_0 = i2_1-1;
				int i2_2 = i2_1 + 1;
				int i3_1 = (int)((edgeAccuracies[typeIndex][charge][0].length-1)*accL);
				int i3_2 = i3_1 + 1;
				if(i3_2 > edgeAccuracies[typeIndex][charge][0].length-1) i3_2 = i3_1-1;
				
				
				float acct0 =0 , acct1 = 0; 

				float acc1 = edgeAccuracies[typeIndex][charge][i1][i2_1][i3_1][0];
				float acc2 = 0, acc0 = 0;
				
				if(i2_0<0){
					acc0 = 0;
				}else
					acc0 = edgeAccuracies[typeIndex][charge][i1][i2_0][i3_1][0];
				
				if(i2_2 >= edgeAccuracies[typeIndex][charge][i1].length){
					acc2 = accR;
				}else
					acc2 = edgeAccuracies[typeIndex][charge][i1][i2_2][i3_1][0];
				//System.out.println(acc1 + "\t" + acc2);
				acct0 = acc1 + (acc2-acc1)/(i2_2-i2_1)*((edgeAccuracies[typeIndex][charge][0].length-1)*accR - i2_1);
				acct0 += acc1 + (acc0-acc1)/(i2_0-i2_1)*((edgeAccuracies[typeIndex][charge][0].length-1)*accR - i2_1);
				acct0/=2;
				
				acc1 = edgeAccuracies[typeIndex][charge][i1][i2_1][i3_2][0];
				acc2 = 0; acc0 = 0;
	
				if(i2_0<0){
					acc0 = 0;
				}else
					acc0 = edgeAccuracies[typeIndex][charge][i1][i2_0][i3_2][0];
				
				if(i2_2 >= edgeAccuracies[typeIndex][charge][i1].length){
					acc2 = accR;
				}else
					acc2 = edgeAccuracies[typeIndex][charge][i1][i2_2][i3_2][0];
				//System.out.println(acc1 + "\t" + acc2);
				acct1 = acc1 + (acc2-acc1)/(i2_2-i2_1)*((edgeAccuracies[typeIndex][charge][0].length-1)*accR - i2_1);
				acct1 += acc1 + (acc0-acc1)/(i2_0-i2_1)*((edgeAccuracies[typeIndex][charge][0].length-1)*accR - i2_1);
				acct1/=2;
				
				acc = acct0 + (acct1-acct0)/(i3_2-i3_1)*((edgeAccuracies[typeIndex][charge][0].length-1)*accL - i3_1);
				
				
			//	acc /= 2;
			 }
			//nacc = r.nullAccuracy;
		}else{
			acc = accR * accL;// nacc = 0;
		}
		
		LR = 0;
		
		acc = Math.min(acc, 1-r.nullAccuracy);// for stability
		acc = Math.max(acc, 0);
		//accuracy = acc;

		
		if(!this.r.isEmptyEnzymaticNode && !this.l.isEmptyEnzymaticNode){
			if(edgeDeviationAccuracies !=null && edgeDeviationAccuracies[typeIndex]!=null && edgeDeviationAccuracies[typeIndex][charge]!=null && edgeDeviationAccuracies[typeIndex][charge][0]!=null  && massDeviationIndex >=0){
				float pp1 = edgeDeviationAccuracies[typeIndex][charge][0][Math.min(edgeDeviationAccuracies[typeIndex][charge][0].length-1, minAANum)][massDeviationIndex];
				float pp2 = edgeDeviationAccuracies[typeIndex][charge][1][Math.min(edgeDeviationAccuracies[typeIndex][charge][0].length-1, minAANum)][massDeviationIndex];
				
				acc = pp1*acc/(pp1*acc + pp2*(1-acc));
				LR += Math.log(pp1/pp2);
			}
		}
		
		if(isPTM()) LR += Math.log(0.05f);
		
		accuracy = Math.max(acc, accL*accR);
		
		float factor = 0.3f;
	//	if(this.isGap) factor /= 3f;
		
		accuracy = Math.min(accuracy, accR + (float)Math.sqrt(accR * (1-accR)) * factor);//
		accuracy = Math.min(accuracy, accL + (float)Math.sqrt(accL * (1-accL)) * factor);// 
		
	}
	
	@Override
	public float getLR(){
		return this.r.getLR() + LR;
	}
	
	public ArrayList<ArrayList<AminoAcid>> getModifiedAAList(){
		return modTable.get(NominalMass.toNominalMass(this.getMass()));
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
		//	if(isCompositeMass) s += "*";
		}
		//System.out.println(aaResidueTable);
		return s.toString();
	}
	
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

}
