package parser;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.Collections;

import msutil.Peptide;

public class MSGFDBParser {
	public static PSMList<PSM> parse(String fileName){
		return parse(fileName, false);
	}
	
	public static PSMList<PSM> parse(String fileName, boolean useFDRScore)
	{
		PSMList<PSM> psmList = new PSMList<PSM>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String s;
		int fdrColumn = -1;
		int pepfdrColumn = -1;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(s.startsWith("#") || s.length() == 0){
				if(useFDRScore){
					for(int i=0; i<token.length;i++){
						String t = token[i];
						if(t.equals("FDR")) fdrColumn = i;
						if(t.equals("PepFDR")) pepfdrColumn = i;
					}
				}
				continue;
				// || !Character.isDigit(s.charAt(0)))
			}
				
			
			//if(token.length < 14)
			//	continue;
			String specFileName = token[0].substring(1+token[0].lastIndexOf(System.getProperty("file.separator")));
			int scanNum = Integer.parseInt(token[1]);
//			String method = token[2];
			float precursorMz = Float.parseFloat(token[4]);
			int charge = Integer.parseInt(token[6]);
			String peptideStr = token[7].substring(token[7].indexOf('.')+1, token[7].lastIndexOf('.'));
			String protein = token[8];
//			int msgfScore = Integer.parseInt(token[7]);
			int peptideScore = Integer.parseInt(token[10]);
			float specProb = Float.parseFloat(token[11]);
			PSM psm = new PSM();
			char precedingResidue = token[7].charAt(0);
			char succeedingResidue = token[7].charAt(token[7].length()-1);
			if(precedingResidue == '.')  precedingResidue = '*';
			if(succeedingResidue == '.')  succeedingResidue = '*';
			
			psm.specFileName(specFileName).scanNum(scanNum).precursorMz(precursorMz).charge(charge).peptide(new Peptide(peptideStr))
			.protein(protein).rawScore(peptideScore).probScore(specProb).precedingResidue(precedingResidue).succeedingResidue(succeedingResidue);
			
			if(useFDRScore && fdrColumn>=0 && pepfdrColumn>=0) {
				psm.score("F", Float.parseFloat(token[fdrColumn]));
				psm.score("P", Float.parseFloat(token[pepfdrColumn]));
			}
			
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
