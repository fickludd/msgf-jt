package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;

import msutil.Peptide;
import msutil.SpectraContainer;
import msutil.SpectraMap;
import msutil.Spectrum;
import msutil.SpectrumAccessorByScanNum;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraMap;

public class AnnotatedSpecGenerator {
	public static void main(String argv[]) throws Exception
	{
		File resultFile = null;
		File specDir = null;
		File outputFile = null;
		int scanNumCol = -1;
		int specFileCol = -1;
		int peptideCol = -1;
		int scoreCol = -1;
		boolean isGreaterBetter = false;
		float threshold;
		if(isGreaterBetter)
			threshold = Float.MIN_VALUE;
		else
			threshold = Float.MAX_VALUE;
		boolean uniquePeptide = false;
		boolean hasHeader = true;
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		int i=0;
		while(i<argv.length)
		{
     		if(argv[i].equalsIgnoreCase("-r"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				resultFile = new File(argv[i+1]);
				if(!resultFile.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist.");
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-d"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				specDir = new File(argv[i+1]);
				if(!specDir.exists() || !specDir.isDirectory())
					printUsageAndExit(argv[i+1] + " doesn't exist or not a directory!");
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-o"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				outputFile = new File(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-n"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				scanNumCol = Integer.parseInt(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-f"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				specFileCol = Integer.parseInt(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-p"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				peptideCol = Integer.parseInt(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-s"))
			{
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				scoreCol = Integer.parseInt(argv[i+1]);
				if(argv[i+2].equalsIgnoreCase("1"))
					isGreaterBetter = true;
				i += 3;
			}
     		else if(argv[i].equalsIgnoreCase("-t"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				threshold = Float.parseFloat(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-u"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				if(argv[i+1].equalsIgnoreCase("1"))
					uniquePeptide = true;
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-h"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				if(argv[i+1].equalsIgnoreCase("0"))
					hasHeader = false;
				i += 2;
			}
//			else if(argv[i].equalsIgnoreCase("-fixMod"))
//			{
//				// 0: No mod, 1: Carbamidomethyl C, 2: Carboxymethyl C
//				if(argv[i+1].equalsIgnoreCase("0"))
//					aaSet = AminoAcidSet.getStandardAminoAcidSet();
//				else if(argv[i+1].equalsIgnoreCase("1"))
//					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//				else if(argv[i+1].equalsIgnoreCase("2"))
//					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarboxymethylatedCys();
//				else
//					printUsageAndExit("Illigal -fixMod parameter: " + argv[i+1]);
//				i += 2;
//			}     		
     		else
     			printUsageAndExit("Illegal parameter!");
		}
		
		if(resultFile == null)
			printUsageAndExit("-r resultFileName is missing");
		if(outputFile == null)
			printUsageAndExit("-o outputFileName is missing");
		if(specDir == null)
			printUsageAndExit("-d specDir is missing");
		if(scanNumCol < 0)
			printUsageAndExit("-n scanNumCol is missing");
		if(specFileCol < 0)
			printUsageAndExit("-f specFileCol is missing");
		if(peptideCol < 0)
			printUsageAndExit("-p peptideCol is missing");
		if(scoreCol < 0)
			printUsageAndExit("-s scoreCol 0/1 is missing");
		
		generateAnnotatedSpectra(resultFile, specDir, outputFile, scanNumCol, specFileCol, peptideCol, scoreCol, isGreaterBetter,
				threshold, uniquePeptide, hasHeader);
		
	}
	
	public static void printUsageAndExit(String message)
	{
		System.out.println(message);
		System.out.println("usage: java AnnotatedSpecGenerator\n" +
				"\t-r resultFileName\n" +
				"\t-d specDir\n" +
				"\t-o outputFileName\n" +
				"\t-f specFileCol\n" +
				"\t-n scanNumCol\n" +
				"\t-p peptideCol\n" +
				"\t-s scoreCol 0/1 (0: smaller is better, 1: larger is better)\n" +
				"\t[-t threshold] \n" +
				"\t[-u 0/1] (0: one spectrum per peptide, 1: no restriction (default))\n" +
				"\t[-h 0/1] (0: no header, 1: header (default))\n"
//				"\t[-fixMod 0/1/2] (0: NoCysteineProtection, 1: CarbamidomethyC (default), 2: CarboxymethylC)\n"
				);
		System.exit(-1);
	}
	
	public static void generateAnnotatedSpectra(File resultFile, File specDir, File outputFile, int scanNumCol, int specFileCol, int pepCol,
			int scoreCol, boolean isGreaterBetter, float threshold, boolean uniquePeptide, boolean hasHeader) throws Exception
	{
		String s;
		Hashtable<String, String> pepTable = null;
		ArrayList<String> resultList = null;
		if(uniquePeptide)
			pepTable = new Hashtable<String, String>();
		else
			resultList = new ArrayList<String>();
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
		if(hasHeader)
			out.println(in.readLine());
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length <= scanNumCol || token.length <= specFileCol || token.length <= pepCol || token.length <= scoreCol)
				continue;
			String pep = token[pepCol];
			if(pep.matches("[A-Z]\\.[A-Z]+\\.[A-Z]"))
				pep = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
			float score = Float.parseFloat(token[scoreCol]);
			
			if(!isGreaterBetter && score < threshold || isGreaterBetter && score > threshold)
			{
				if(uniquePeptide)
				{
					if(pepTable.get(pep) == null)
						pepTable.put(pep, s);
					else
					{
						String[] token2 = pepTable.get(pep).split("\t");
						float existingScore = Float.parseFloat(token2[scoreCol]);
						if(score < existingScore)
							pepTable.put(pep, s);
					}
				}
				else
					resultList.add(s);
			}
		}
		
		Iterator<String> itr;
		if(uniquePeptide)
			itr = pepTable.values().iterator();
		else
			itr = resultList.iterator();

		HashMap<String,SpectrumAccessorByScanNum> specAccessorMap = new HashMap<String,SpectrumAccessorByScanNum>(); 
		SpectraContainer container = new SpectraContainer();
		while(itr.hasNext())
		{
			String str = itr.next();
			String[] token = str.split("\t");
			
			String specFileName = token[specFileCol];
			SpectrumAccessorByScanNum specMap = specAccessorMap.get(specFileName);
			if(specMap == null)
			{
				String ext = specFileName.substring(specFileName.lastIndexOf('.'));
				if(ext.equalsIgnoreCase(".mzXML"))
					specMap = new MzXMLSpectraMap(specDir.getPath()+File.separator+specFileName);
				else if(ext.equalsIgnoreCase(".mgf"))
					specMap = new SpectraMap(specDir.getPath()+File.separator+specFileName, new MgfSpectrumParser());
				else
				{
					System.out.println("Unrecognized spectrum format: " + specFileName);
					System.exit(-1);
				}
				specAccessorMap.put(specFileName, specMap);
			}
			
			String pep = token[pepCol];
			if(pep.matches(".*\\.[A-Z]+\\..*"))
				pep = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
			int scanNum = Integer.parseInt(token[scanNumCol]);
			
			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
			if(spec == null)
			{
				System.out.println(specFileName+":"+scanNum+" is not available!");
				System.exit(-1);
			}
			spec.setAnnotation(new Peptide(pep));
			container.add(spec);
		}
		container.outputMgfFile(outputFile.getPath());
		System.out.println("Done");		
	}
	
}
