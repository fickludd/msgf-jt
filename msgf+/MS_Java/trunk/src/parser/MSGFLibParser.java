package parser;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.Collections;

import msutil.Peptide;

public class MSGFLibParser {
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
			if(s.startsWith("#") || s.length() == 0)// || !Character.isDigit(s.charAt(0)))
				continue;
			String[] token = s.split("\t");
			//		ps.println("#SpecFileName\tScan#\tCharge\tPrecursorMz\tPeptide\tProtein\tpValue");

			String specFileName = token[0].substring(1+token[0].lastIndexOf(System.getProperty("file.separator")));
			int scanNum = Integer.parseInt(token[1]);
//			String method = token[2];
			float precursorMz = Float.parseFloat(token[3]);
			int charge = Integer.parseInt(token[2]);
			String peptideStr = token[4];//.substring(token[4].indexOf('.')+1, token[4].lastIndexOf('.'));
			String protein = token[5];
//			int msgfScore = Integer.parseInt(token[7]);
		//	int peptideScore = Integer.parseInt(token[10]);
			float specProb = Float.parseFloat(token[8]);
			
			
			PSM psm = new PSM();
			psm.specFileName(specFileName).scanNum(scanNum).precursorMz(precursorMz).charge(charge).peptide(new Peptide(peptideStr))
			.protein(protein).probScore(specProb).score("M", Float.parseFloat(token[6])).score("S", Float.parseFloat(token[7])).score("IsEPvalue", Float.parseFloat(token[9]));
			psmList.add(psm);
		}
		return psmList;
	}
	
	public static PSMList<PSM> getMergedResults(String dirName, final String suffix)
	{
		File dir = new File(dirName);
		if(!dir.isDirectory())
			return null;
		class SuffixFileFilter implements FileFilter {
			public boolean accept(final File pathname) {
				if(pathname.getName().endsWith(suffix))
					return true;
				else
					return false;
			}
		}
		PSMList<PSM> psmList = new PSMList<PSM>();
		for(File f : dir.listFiles(new SuffixFileFilter()))
		{
			psmList.addAll(parse(f.getPath()));
		}
		Collections.sort(psmList, new PSM.PSMSpecNumComparator());
		return psmList;
	}
}
