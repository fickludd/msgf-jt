package parser;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.Collections;

import msutil.Peptide;

public class SLGFParser {
	
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
			

			if(s.length() == 0) continue;
			
			String[] token = s.split("\t");
			
			if(s.startsWith("Scan#")){
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
				//continue;
			//"Scan#\t#SpectrumFile\tAnnotation\tOrigAnnotation\tProtein\tdbIndex\tnumMods\tmatchOrientation\tstartMass\tCharge\tMQScore\tp-value\tisDecoy\tStrictEnvelopeScore\tUnstrictEvelopeScore\tMZ\tFDR";

			String specFileName = token[1].substring(1+token[1].lastIndexOf(System.getProperty("file.separator")));;
			//System.out.println(s);
			int scanNum = Integer.parseInt(token[0]);
//			String method = token[2];
			float precursorMz = Float.parseFloat(token[15]);
			int charge = Integer.parseInt(token[9]);
			String peptideStr = token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.'));
			
			float pval = Float.parseFloat(token[10]);
			
			PSM psm = new PSM();
			psm.specFileName(specFileName).scanNum(scanNum).charge(charge).peptide(new Peptide(peptideStr))
			.probScore(pval).precursorMz(precursorMz).precedingResidue('*').succeedingResidue('*');
			
			if(useFDRScore && fdrColumn>=0 && pepfdrColumn>=0) {
				psm.score("F", Float.parseFloat(token[fdrColumn]));
				psm.score("P", Float.parseFloat(token[pepfdrColumn]));
			}
			
			if(token[12].equals("1")) psm.score("DECOY", 1);
			else psm.score("DECOY", 0);
			
			//if(token.length > 7)
			//	psm.rawScore(Float.parseFloat(token[7]));
			
			
			
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
