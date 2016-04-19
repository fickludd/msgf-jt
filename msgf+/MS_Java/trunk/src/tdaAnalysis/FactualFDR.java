package tdaAnalysis;

import msgf.Tolerance;
import msutil.Composition;
import parser.PSM;
import parser.PSMList;

public class FactualFDR {
	
	private boolean proteinFIlter = false;
	private boolean pmFilter = false;
	private boolean specFilter = false;
	private boolean isISB, isMSGFDB;
	
	private String[] proteinKeys = null;
	private int maxSN = 0;
	private String specKey = null;
	private Tolerance pmTolerance = null;
	
	private PSMList<PSM> positivePSMs = null;
	private PSMList<PSM> falsePSMs = null;
	private PSMList<PSM> passedPSMs = null;
	
	public FactualFDR(PSMList<PSM> positivePSMs, boolean isISB, boolean isMSGFDB){
		this.positivePSMs = positivePSMs.getDistinctiveSpectralSet();	
		this.isISB = isISB;
		this.isMSGFDB = isMSGFDB;
	}
	
	public FactualFDR(){
		
	}
	
	public FactualFDR setProteinFilter(String[] keys){
		proteinFIlter = true;
		proteinKeys = keys;
		return this;
	}
	
	public FactualFDR setSpecFilter(int maxSN){ // inclusive
		specFilter = true;
		this.maxSN = maxSN;
		return this;
	}
	
	public FactualFDR setSpecFilter(String key){ // inclusive
		specFilter = true;
		this.specKey = key;
		return this;
	}
	
	public FactualFDR setPMFilter(Tolerance tol){
		pmFilter = true;
		pmTolerance = tol;
		return this;
	}
	
	public void calculate(){
		falsePSMs = new PSMList<PSM>();
		passedPSMs = new PSMList<PSM>();
		for(PSM psm : positivePSMs){
			boolean pass = true;
			
			if(pass && proteinFIlter){
				boolean proPass = false;
				
				for(String proteinKey : proteinKeys){
					if(psm.getProtein().startsWith(proteinKey)){
						proPass = true;
						break;
					}
				}
				
				pass = proPass;
				//if(!pass) System.out.println(psm.getProtein());
			}
			if(pass && specFilter){
				if(psm.getTitle() != null && specKey !=null && !psm.getTitle().contains(specKey)){
					pass = false;
				}
				
				if(specKey == null && psm.getScanNum() > maxSN){
					pass = false;
				}
			}
			if(pass && pmFilter){
				float pm = (psm.getPrecursorMz() - (float)Composition.PROTON) * psm.getCharge(); 
				boolean pmPass = false;
				if(Math.abs(pm - psm.getPeptide().getParentMass()) <= pmTolerance.getToleranceAsDa(pm)){
					pmPass = true;
				//	System.out.println(pm + "\t" + psm.getPeptide().getParentMass());
				}
			if(!pmPass){
					float pmt = (float) (pm - Composition.ISOTOPE);
					if(Math.abs(pmt - psm.getPeptide().getParentMass()) <= pmTolerance.getToleranceAsDa(pmt)){
						pmPass = true;
					}
				}
				if(!pmPass){
					float pmt = (float) (pm - 2*Composition.ISOTOPE);
					if(Math.abs(pmt - psm.getPeptide().getParentMass()) <= pmTolerance.getToleranceAsDa(pmt)){
						pmPass = true;
					}
				}
				if(!pmPass){
					float pmt = (float) (pm - Composition.ISOTOPE2);
					if(Math.abs(pmt - psm.getPeptide().getParentMass()) <= pmTolerance.getToleranceAsDa(pmt)){
						pmPass = true;
					}
				}
				
				
				if(!pmPass){
					float pmt = (float) (pm + Composition.ISOTOPE);
					if(Math.abs(pmt - psm.getPeptide().getParentMass()) <= pmTolerance.getToleranceAsDa(pmt)){
						pmPass = true;
					}
				}
				if(!pmPass){
					float pmt = (float) (pm + 2*Composition.ISOTOPE);
					if(Math.abs(pmt - psm.getPeptide().getParentMass()) <= pmTolerance.getToleranceAsDa(pmt)){
						pmPass = true;
					}
				}
				if(!pmPass){
					float pmt = (float) (pm + Composition.ISOTOPE2);
					if(Math.abs(pmt - psm.getPeptide().getParentMass()) <= pmTolerance.getToleranceAsDa(pmt)){
						pmPass = true;
					}
				}
				
				if(!pmPass){
					if(Math.abs(pm - psm.getPeptide().getParentMass()) > 10){ // weird XTandem
						pmPass = true;
					}
				}
				
				
				//if(!pmPass) System.out.println(psm.getScanNum() + "\t" + psm.getPeptideStr() + "\t" + pm + "\t" + psm.getPeptide().getParentMass());
				
				pass = pmPass;
				
				
			}
			if(!pass){
				//System.out.println(psm+"\t"+(Math.abs((psm.getPrecursorMz() - (float)Composition.PROTON) * psm.getCharge() - psm.getPeptide().getParentMass())) + "\t" +  pmTolerance.getToleranceAsDa((psm.getPrecursorMz() - (float)Composition.PROTON) * psm.getCharge()));
				falsePSMs.add(psm);
			}else passedPSMs.add(psm);
		}
	}
	
	public PSMList<PSM> getPassedPSMs() {return passedPSMs;}
	
	public PSMList<PSM> getFalsePSMs() {return falsePSMs;}
	
	public int getNPositivePeps(){
		return positivePSMs.getDistinctivePeptideSet().size();
	}
	
	public int getNfalsePeps(){
		return falsePSMs.getDistinctivePeptideSet().size();
	}
	
	public int getNPositive(){
		return positivePSMs.getDistinctiveSpectralSet().size();
	}
	
	public int getNfalse(){
		return falsePSMs.getDistinctiveSpectralSet().size();
	}

	
}
