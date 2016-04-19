package parser;

import java.io.FileNotFoundException;

import msutil.Composition;
import msutil.Peptide;

public class XTandemParser {
	public static PSMList<PSM> parse(String fileName)
	{
		PSMList<PSM> psmList = new PSMList<PSM>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String s;
		int chargeFieldOffset = 0;
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpectrumFile") || s.length() == 0){
				if(s.contains("Charge")) chargeFieldOffset = 1;
				continue;
				// || !Character.isDigit(s.charAt(0)))
			}
				
			String[] token = s.split("\t");
			
			String specFileName = token[0];
			int scanNum = Integer.parseInt(token[1]);
//			String method = token[2];
			
			String peptideStr = token[2];//token[7].substring(token[7].indexOf('.')+1, token[7].lastIndexOf('.'));
			String protein = token[3];
			int charge = 1;// TODO
			if(chargeFieldOffset > 0) charge = Integer.parseInt(token[4]);
			
			float precursorMz = (Float.parseFloat(token[4+chargeFieldOffset]));
			precursorMz = (float) ((precursorMz - Composition.PROTON + charge * Composition.PROTON)/charge);
			
			
			
//			int msgfScore = Integer.parseInt(token[7]);
			float E_value = Float.parseFloat(token[7+chargeFieldOffset]);
			
			float peptideScore = Float.parseFloat(token[8+chargeFieldOffset]);
			PSM psm = new PSM();
			
			
			
			psm.specFileName(specFileName).scanNum(scanNum).precursorMz(precursorMz).peptide(new Peptide(peptideStr))
			.protein(protein).rawScore(peptideScore).probScore(E_value).charge(charge);
			
			float pm = (psm.getPrecursorMz() - (float)Composition.PROTON) * psm.getCharge(); 
			if(Math.abs(pm - psm.getPeptide().getParentMass()) < 3) // TODO remove!
				
				psmList.add(psm);
		}
		return psmList;
	}
}
