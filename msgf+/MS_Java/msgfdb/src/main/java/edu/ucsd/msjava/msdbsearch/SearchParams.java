package edu.ucsd.msjava.msdbsearch;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.DBSearchIOFiles;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpecFileFormat;
import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.IntRangeParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.params.ToleranceParameter;

public class SearchParams {
	private List<DBSearchIOFiles> dbSearchIOList;
	private File databaseFile;
	private Tolerance leftParentMassTolerance; 
	private Tolerance rightParentMassTolerance; 
	private int minIsotopeError;
	private int maxIsotopeError;
	private Enzyme enzyme;
	private int numTolerableTermini;
	private ActivationMethod activationMethod;
	private InstrumentType instType;
	private Protocol protocol;
	private AminoAcidSet aaSet;
	private int numMatchesPerSpec;
	private int startSpecIndex;
	private int endSpecIndex;
	private boolean useTDA;
//	private boolean showFDR;
//	private boolean showDecoy;
	private int minPeptideLength;
	private int maxPeptideLength;
	private int maxNumVariatsPerPeptide;
	private int minCharge;
	private int maxCharge;
	private double intensityThreshold;
	private int numThreads;
	private boolean replicateMergedResults;
	private boolean doNotDseEdgeScore;
	private File dbIndexDir;
	private boolean outputAdditionalFeatures;
	private int minNumPeaksPerSpectrum;
	
	
	public SearchParams()	{}
	
	public List<DBSearchIOFiles> getDBSearchIOList() {
		return dbSearchIOList;
	}

	public File getDatabaseFile() {
		return databaseFile;
	}

	public Tolerance getLeftParentMassTolerance() {
		return leftParentMassTolerance;
	}

	public Tolerance getRightParentMassTolerance() {
		return rightParentMassTolerance;
	}

	public int getMinIsotopeError() {
		return minIsotopeError;
	}

	public int getMaxIsotopeError() {
		return maxIsotopeError;
	}
	
	public Enzyme getEnzyme() {
		return enzyme;
	}

	public int getNumTolerableTermini() {
		return numTolerableTermini;
	}

	public ActivationMethod getActivationMethod() {
		return activationMethod;
	}

	public InstrumentType getInstType() {
		return instType;
	}

	public Protocol getProtocol() {
		return protocol;
	}

	public AminoAcidSet getAASet() {
		return aaSet;
	}

	public int getNumMatchesPerSpec() {
		return numMatchesPerSpec;
	}

	public int getStartSpecIndex() {
		return startSpecIndex;
	}

	public int getEndSpecIndex() {
		return endSpecIndex;
	}

	public boolean useTDA() {
		return useTDA;
	}

//	public boolean showFDR() {
//		return showFDR;
//	}
//
//	public boolean showDecoy() {
//		return showDecoy;
//	}

	public int getMinPeptideLength() {
		return minPeptideLength;
	}

	public int getMaxPeptideLength() {
		return maxPeptideLength;
	}

	public int getMaxNumVariatsPerPeptide() {
		return maxNumVariatsPerPeptide;
	}
	
	public int getMinCharge() {
		return minCharge;
	}

	public int getMaxCharge() {
		return maxCharge;
	}

	public double getIntensityThreshold() {
		return intensityThreshold;
	}

	public int getNumThreads() {
		return numThreads;
	}

	public boolean replicateMergedResults() {
		return replicateMergedResults;
	}

	public boolean doNotDseEdgeScore() {
		return doNotDseEdgeScore;
	}

	public File getDBIndexDir() {
		return dbIndexDir;
	}
	
	public boolean outputAdditionalFeatures()
	{
		return outputAdditionalFeatures;
	}
	
	public int getMinNumPeaksPerSpectrum()
	{
		return minNumPeaksPerSpectrum;
	}
	
	public String parse(ParamManager paramManager)
	{
		// Spectrum file
		FileParameter specParam = paramManager.getSpecFileParam();
		File specPath = specParam.getFile();
		
		dbSearchIOList = new ArrayList<DBSearchIOFiles>();
		
		if(!specPath.isDirectory())
		{
			// Spectrum format
			SpecFileFormat specFormat = (SpecFileFormat)specParam.getFileFormat();
			// Output file
			File outputFile = paramManager.getOutputFileParam().getFile();
			if(outputFile == null)
			{
				String outputFilePath = specPath.getPath().substring(0, specPath.getPath().lastIndexOf('.'))+".mzid";
				outputFile = new File(outputFilePath);
//				if(outputFile.exists())
//					return outputFile.getPath() + " already exists!";
			}
			
			dbSearchIOList = new ArrayList<DBSearchIOFiles>();
			dbSearchIOList.add(new DBSearchIOFiles(specPath, specFormat, outputFile));
		}
		else	// spectrum directory
		{
			dbSearchIOList = new ArrayList<DBSearchIOFiles>();
			for(File f : specPath.listFiles())
			{
				SpecFileFormat specFormat = SpecFileFormat.getSpecFileFormat(f.getName());
				if(specParam.isSupported(specFormat))
				{
					String outputFileName = f.getName().substring(0, f.getName().lastIndexOf('.'))+".mzid";
					File outputFile = new File(outputFileName);
//					if(outputFile.exists())
//						return outputFile.getPath() + " already exists!";
					dbSearchIOList.add(new DBSearchIOFiles(f, specFormat, outputFile));
				}
			}
		}
		
		// DB file
		databaseFile = paramManager.getDBFileParam().getFile();
		
		// PM tolerance
		ToleranceParameter tol = ((ToleranceParameter)paramManager.getParameter("t"));
		leftParentMassTolerance = tol.getLeftTolerance();
		rightParentMassTolerance = tol.getRightTolerance();
		
		int toleranceUnit = paramManager.getIntValue("u");
		if(toleranceUnit != 2)
		{
			boolean isTolerancePPM;
			if(toleranceUnit == 0)
				isTolerancePPM = false;
			else 
				isTolerancePPM = true;
			leftParentMassTolerance = new Tolerance(leftParentMassTolerance.getValue(), isTolerancePPM);
			rightParentMassTolerance = new Tolerance(rightParentMassTolerance.getValue(), isTolerancePPM);
		}
		
		IntRangeParameter isotopeParam = (IntRangeParameter)paramManager.getParameter("ti");
		this.minIsotopeError = isotopeParam.getMin();
		this.maxIsotopeError = isotopeParam.getMax();
		
		if(rightParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f || leftParentMassTolerance.getToleranceAsDa(1000, 2) >= 0.5f)
			minIsotopeError = maxIsotopeError = 0;
		
		enzyme = paramManager.getEnzyme();
		numTolerableTermini = paramManager.getIntValue("ntt");
		activationMethod = paramManager.getActivationMethod();
		instType = paramManager.getInstType();
		if(activationMethod == ActivationMethod.HCD && instType != InstrumentType.HIGH_RESOLUTION_LTQ && instType != InstrumentType.QEXACTIVE)
			instType = InstrumentType.QEXACTIVE;	// by default use Q-Exactive model for HCD
		
		protocol = paramManager.getProtocol();
		
		aaSet = null;
		File modFile = paramManager.getModFileParam().getFile();
		if(modFile == null)
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		else
		{
			String modFileName = modFile.getName();
			String ext = modFileName.substring(modFileName.lastIndexOf('.')+1);
			if(ext.equalsIgnoreCase("xml"))
				aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(modFile.getPath());
			else
				aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
			if(protocol == Protocol.NOPROTOCOL)
			{
				if(aaSet.containsITRAQ())
				{
					if(aaSet.containsPhosphorylation())
						protocol = Protocol.ITRAQPHOSPHO;
					else
						protocol = Protocol.ITRAQ;
				}
				else if(aaSet.containsTMT())
				{
					protocol = Protocol.TMT;
				}
				else
				{
					if(aaSet.containsPhosphorylation())
						protocol = Protocol.PHOSPHORYLATION;
					else
						protocol = Protocol.NOPROTOCOL;
				}
			}
			
		}
		
		numMatchesPerSpec = paramManager.getIntValue("n");
		
		startSpecIndex = ((IntRangeParameter)paramManager.getParameter("index")).getMin();
		endSpecIndex = ((IntRangeParameter)paramManager.getParameter("index")).getMax();
		
		useTDA = paramManager.getIntValue("tda") == 1 ? true : false;
//		showFDR = paramManager.getIntValue("showQValue") == 1 ? true : false;
//		showDecoy = paramManager.getIntValue("showDecoy") == 1 ? true : false;
		outputAdditionalFeatures = paramManager.getIntValue("addFeatures") == 1 ? true : false;
		
		minPeptideLength = paramManager.getIntValue("minLength");
		maxPeptideLength = paramManager.getIntValue("maxLength");
		
		maxNumVariatsPerPeptide = paramManager.getIntValue("iso");
		
		if(minPeptideLength > maxPeptideLength)
		{
			return "MinPepLength must not be larger than MaxPepLength";
		}
		
		minCharge = paramManager.getIntValue("minCharge");
		maxCharge = paramManager.getIntValue("maxCharge");
		if(minCharge > maxCharge)
		{
			return "MinCharge must not be larger than MaxCharge";
		}
		
		intensityThreshold = paramManager.getFloatValue("intensityThreshold");
		
		numThreads = paramManager.getIntValue("thread");
		doNotDseEdgeScore = paramManager.getIntValue("edgeScore") == 1 ? true : false;
		
		dbIndexDir = paramManager.getFile("dd");
		
		minNumPeaksPerSpectrum = paramManager.getIntValue("minNumPeaks");
		
		return null;
	}
	
	@Override
	public String toString()
	{
		StringBuffer buf = new StringBuffer();

//		buf.append("Spectrum File(s):\n");
//		for(DBSearchIOFiles ioFile : this.dbSearchIOList)
//		{
//			buf.append("\t"+ioFile.getSpecFile().getAbsolutePath()+"\n");
//		}
//		buf.append("Database File: " + this.databaseFile.getAbsolutePath() + "\n");

		buf.append("\tPrecursorMassTolerance: ");
		if(leftParentMassTolerance.equals(rightParentMassTolerance))
			buf.append(leftParentMassTolerance);
		else
			buf.append("["+leftParentMassTolerance+","+rightParentMassTolerance+"]");
		buf.append("\n");
		
		buf.append("\tIsotopeError: " + this.minIsotopeError + "," + this.maxIsotopeError + "\n");
		buf.append("\tTargetDecoyAnalysis: " + this.useTDA + "\n");
		buf.append("\tFragmentationMethod: " + this.activationMethod + "\n");
		buf.append("\tInstrument: " + (instType == null ? "null" : this.instType.getName()) + "\n");
		buf.append("\tEnzyme: " + (enzyme == null ? "null" : this.enzyme.getName()) + "\n");
		buf.append("\tProtocol: " + (protocol == null ? "null" : this.protocol.getName()) + "\n");
		buf.append("\tNumTolerableTermini: " + this.numTolerableTermini + "\n");
		buf.append("\tMinPeptideLength: " + this.minPeptideLength + "\n");
		buf.append("\tMaxPeptideLength: " + this.maxPeptideLength + "\n");
		buf.append("\tNumMatchesPerSpec: " + this.numMatchesPerSpec);
		
		return buf.toString();
	}
}