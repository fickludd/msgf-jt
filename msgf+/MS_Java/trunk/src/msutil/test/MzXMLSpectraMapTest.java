package msutil.test;

import msutil.Spectrum;
import parser.MzXMLSpectraMap;

/**
 * Test class for the object of this name
 * @author jung
 *
 */
public class MzXMLSpectraMapTest {

  public static void main(String argv[])
  {
	  test1();
  }
  
  public static void test1()
  {
	    //String fileName = System.getProperty("user.home")+"/Research/Data/HEK293/H293-total-try-a-200ug-2D34-081905-LTQ1-01.mzXML";
	    
	    // Big file
	    String fileName = System.getProperty("user.home")+"/Data/Cyclic/Javier/cyclicPeptides/MTVVI-20_1.medium.mzXML";
	  
	    MzXMLSpectraMap map = new MzXMLSpectraMap(fileName);
	    System.out.println("Total number of scans: " + map.getScanCount());
	    
	    for (int i = 1; i <= map.getScanCount(); i++) {
	      Spectrum s = map.getSpectrumByScanNum(i);
	      if (s != null) {
	        System.out.println("Scan number found " + s.getScanNum());
	      }
	      else {
	        System.out.println("Scan number not found " + i);
	      }
	    }	  
  }
}