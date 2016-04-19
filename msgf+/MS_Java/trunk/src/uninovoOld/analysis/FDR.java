package uninovoOld.analysis;

import parser.MODAParser;
import parser.PSM;
import parser.PSMList;

public class FDR {
	public static void main(String[] args){
		PSMList<PSM> PSMs = new PSMList<PSM>();//MSGFDBParser.parse("/home/kwj/UniNovo/MSGFDB.txt");
		//PSMs.addAll(UniNovoParser.parse("/home/kwj/UniNovo/UniCIDNoM16.txtmodMSGF"));
		//PSMs.addAll(MSGFDBParser.parse("/home/kwj/UniNovo/MSGFDBNoM16REV.txt"));
		//PSMs.addAll(MSGFDBParser.parse("/home/kwj/UniNovo/MSGFDBREV.txt"));
		//PSMs.addAll(UniNovoParser.parse("/home/kwj/UniNovo/UniCID.txtmodMSGF"));
		//PSMs.addAll(UniNovoParser.parse("/home/kwj/UniNovo/UniForMSGDCID.txtmodMSGF"));
		
		PSMs.addAll(MODAParser.parse("/home/kwj/UniNovo/MODA-3aba564e-group_by_spectrum-main.tsv"));
		
		
		PSMs = PSMs.getDistinctiveSpectralSet();
		
		//System.out.println(PSMs.size());
		PSMList<PSM> targetPSMList = new PSMList<PSM>();
		PSMList<PSM> decoyPSMList = new PSMList<PSM>();
		
		for(PSM psm : PSMs){
			if(psm.getProtein().contains("REV_")) decoyPSMList.add(psm);
			else targetPSMList.add(psm);
		}
		
		PSMList<PSM> qPSMs = PSMList.selectUsingFDR(targetPSMList, decoyPSMList, 0.01f);
		qPSMs = qPSMs.getDistinctivePeptideSet();
		int p = 0;
		for(PSM psm : qPSMs){
			//System.out.println(new Peptide("K.EGGANLLTTVVNKGPLGQITTGKGTTKQVVIVNPK.T"));
			if(psm.getPeptide().isModified()) p++;
			//if(psm.getPeptideStr().contains("m")) p++;
		//	if(psm.getPeptideStr().contains("s")) p++;
		//	if(psm.getPeptideStr().contains("t")) p++;
		//	if(psm.getPeptideStr().contains("y")) p++;
		}
		
		System.out.println(qPSMs.size() + "\t" + p);
		
	}
}
