package parser;

import java.io.FileNotFoundException;

import msutil.Peptide;

public class UniNovoParser {
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
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#SpectrumFile") || s.length() == 0){
				continue;
			}
				
			String[] token = s.split("\t");
			
			String specFileName = token[0];
			int scanNum = Integer.parseInt(token[1]);
//			String method = token[2];
			
			String peptideStr = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'));
			String protein = token[8];
			int charge = Integer.parseInt(token[4]);
			
			float precursorMz = (Float.parseFloat(token[3]));
			
			
//			int msgfScore = Integer.parseInt(token[7]);
			float p_value = Float.parseFloat(token[token.length-1]);
			
			PSM psm = new PSM();
			
			
			
			psm.specFileName(specFileName).scanNum(scanNum).precursorMz(precursorMz).peptide(new Peptide(peptideStr))
			.protein(protein).probScore(p_value).charge(charge);
			
			psmList.add(psm);
		}
		return psmList;
	}
}
