package parser;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.Collections;

import msutil.Peptide;

public class SWATHParser {
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
			//if(token.length < 14)
				//continue;
			int scanNum = Integer.parseInt(token[0].split(": ")[1]);
//			String method = token[2];
			String peptideStr = token[3];
			//System.out.println(s);
			int charge = Integer.parseInt(token[5]);
			//System.out.println(charge);
			float specProb = 1 -  Float.parseFloat(token[7]);

			
			PSM psm = new PSM();
			psm.scanNum(scanNum).peptide(new Peptide(peptideStr)).charge(charge)
			.probScore(specProb);
		
			psm.score("LSN", Integer.parseInt(token[8].split(":")[1].trim()));
			
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

