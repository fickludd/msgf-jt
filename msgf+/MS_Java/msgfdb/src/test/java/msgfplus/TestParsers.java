package msgfplus;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

import org.junit.Test;

import edu.ucsd.msjava.misc.ConvertToMgf;
import edu.ucsd.msjava.misc.MS2ToMgf;
import edu.ucsd.msjava.misc.ParamToTxt;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.PNNLSpectraIterator;
import edu.ucsd.msjava.parser.PNNLSpectrumParser;

public class TestParsers {
	@Test
	public void testReadingIPRG2014Mgf()
	{
		File mgfFile = new File("D:\\Research\\Data\\IPRG2014\\MGF\\C~~data~iPRG 2014~130716_iPRG14_004.raw.-1.mgf");
		SpectraAccessor specAcc = new SpectraAccessor(mgfFile);
		Iterator<Spectrum> itr = specAcc.getSpecItr();
		int numSpecs = 0;
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			numSpecs++;
		}
		System.out.println("NumSpcs: " + numSpecs);
	}
	
	@Test
	public void testReadingTripleTOFFile()
	{
		File mzMLFile = new File("/Users/kims336/Research/Data/ImmunoPeptidomics/MHC_class1_TripleTOF_5600/mzML/carone_L130326_003_pMHC_1-5ug.mzML");
		File mgfFile = new File("/Users/kims336/Research/Data/ImmunoPeptidomics/MHC_class1_TripleTOF_5600/mzML/carone_L130326_003_pMHC_1-5ug.mgf");
		try {
			ConvertToMgf.convert(mzMLFile, mgfFile, false, null, null, -1, -1, -1, false); 
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	@Test
	public void testParamToTxt()
	{
//		File paramFile = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/params/HCD_HighRes_Tryp_iTRAQ.param");
//		ParamToTxt.convert(paramFile);
	}
	
	@Test
	public void ms2ToMgfTest()
	{
//		File ms2File = new File(System.getProperty("user.home")+"/Research/Data/Viktor/103111-Yeast-2hr-01.ANNOTATED.ms2");
//		File mgfFile = new File(System.getProperty("user.home")+"/Research/Data/Viktor/103111-Yeast-2hr-01.ANNOTATED.mgf");
//		try {
//			MS2ToMgf.convert(ms2File, mgfFile);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
	}
		
//	@Test
//	public void testReadingDtaTxt()
//	{
//		int numSpecs = 0;
//		File mgfFile = new File("C:\\cygwin\\home\\kims336\\Data\\QCShewQE\\QC_Shew_13_02_2500ng_B_18Mar13_Jaguar_13-03-11_dta.txt");
//		try {
//			SpectraIterator itr = new PNNLSpectraIterator(mgfFile.getPath());
//			while(itr.hasNext())
//			{
//				Spectrum spec = itr.next();
//				boolean isHighPrecision = spec.isHighPrecision();
//				if(!isHighPrecision)
//				{
//					System.out.println("ScanNum " + spec.getScanNum());
//				}
//				++numSpecs;
//			}
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//		System.out.println("NumSpectra: " + numSpecs);
//	}
	
}
