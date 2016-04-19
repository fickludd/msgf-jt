package msutil;


import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import msgf.Tolerance;

/**
 * Representation of a mass spectrum object.
 * @author Sangtae Kim
 *
 */
public class Spectrum extends ArrayList<Peak> implements Comparable<Spectrum> {

  //this is recommended for Serializable objects
  static final private long serialVersionUID = 1L;

  // required members
  private Peak precursor = null;

  // optional members
  private int startScanNum = -1;
  private int endScanNum = -1;
  private String title = null;
  private Peptide annotation = null;
  private ArrayList<String> seqList = null;	// SEQ fields of mgf spectrum
  private float rt = -1;                    // retention time
  private ActivationMethod activationMethod = null;	// fragmentation method
  private int msLevel = 2;	// ms level
  
/***** CONSTRUCTORS *****/
  /**
   * Empty constructor.
   */
  public Spectrum() {
  }
  
  /**
   * Constructor from a precursor peak.
   * @param precursorPeak the precursor peak.
   */
  public Spectrum(Peak precursorPeak) {
    this.precursor = precursorPeak;
  }
  
  /**
   * Constructor.
   * @param precursorMz              m/z value read from the file.
   * @param charge                     charge read from the file.
   * @param precursorIntensity         intensity of the precursor peak.
   */
  public Spectrum(float precursorMz, int charge, float precursorIntensity) {
    this.precursor = new Peak(precursorMz,  precursorIntensity, charge);
  }
  
  
  
/***** GETTERS ****/
  /**
   * Gets the Peptide object of the annotation.
   * @return the Peptide object representing the annotation or null if not annotated.
   */
  public Peptide getAnnotation()	{ return annotation; }
  
  /**
   * Gets the annotation as a String object.
   * @return the String object of the annotation. null if there is no annotation.
   */
  public String getAnnotationStr()	{ if(annotation != null) return annotation.toStringWithModification(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys()); return null; }
  /**
   * Gets the list of seq fi as String objects.
   * @return the ArrayList of String objeccts of the annotations. null if there is no annotation.
   */
  public ArrayList<String> getSeqList() { return seqList; }
  
  /**
   * Gets the charge of this spectrum.
   * @return the integer charge.
   */
  public int getCharge() { return precursor.getCharge(); }
  
  /**
   * Gets the end scan number of this spectrum.
   * @return the end scan for this spectrum if present. -1 otherwise.
   */
  public int getEndScanNum() { return endScanNum; }
  
  /**
   * Gets the absolute parent mass of this spectrum.
   * @return the mass in Daltons.
   */
  public float getParentMass() { return precursor.getMass(); }

  /**
   * Gets the peptide mass of this spectrum: parentMass-mass(H2O)
   * @return the peptide mass in Daltons.
   */
  public float getPeptideMass() { return precursor.getMass()-(float)Composition.H2O; }
  
  /**
   * Gets the precursor peak member.
   * @return the peak.
   */
  public Peak getPrecursorPeak() { return precursor; }

  /**
   * Gets the scan number of this spectrum.
   * @return the start scan number of this spectrum if present. -1 otherwise.
   */
  public int getScanNum() { return getStartScanNum(); }
  
  /**
   * Gets the start scan number of the spectrum.
   * @return the start scan for this spectrum if present. -1 otherwise.
   */
  public int getStartScanNum() { return startScanNum; }
  
  /**
   * Gets the title of this spectrum.
   * @return the title of this spectrum or null if it was not specified.
   */
  public String getTitle() { return title; }

  /**
   * Returns the retention time in seconds.
   * @return the retention time if set, negative value otherwise.
   */
  public float getRt() { return this.rt; }
  
  /**
   * Returns the activation method.
   * @return the activation method or null if not specified.
   */
  public ActivationMethod getActivationMethod() { return this.activationMethod; }
  
  /**
   * Returns the ms level.
   * @return the ms level.
   */
  public int getMSLevel()	{ return this.msLevel; }
  
/***** SETTERS *****/
  
  /**
   * Sets the annotation with a Peptide object.
   * @param annotation	annotation object
   */
  public void setAnnotation(Peptide annotation) { this.annotation = annotation; }

  /**
   * Sets the annotation with a Sequence object.
   * @param annotation                 annotation object.
   */
  public void addSEQ(String seq) 
  { 
	  if(seqList == null)
		  seqList = new ArrayList<String>();
	  this.seqList.add(seq);
  }

  /**
   * Set the precursor mass.
   * @param precursor 
   */
  public void setPrecursor(Peak precursor) { this.precursor = precursor; }
      
  /**
   * Sets the starting scan number of this spectrum.
   * @param startScanNum         starting scan number from the file.
   */
  public void setStartScanNum(int startScanNum) { this.startScanNum = startScanNum; }

  /**
   * Sets the ending scan number of this spectrum.
   * @param endScanNum         ending scan number from the file.
   */
  public void setEndScanNum(int endScanNum) { this.endScanNum = endScanNum; }
  
  /**
   * Sets the scan number of this spectrum.
   * @param scanNum                    scan number from the file.
   */
  public void setScanNum(int scanNum) { this.startScanNum = scanNum; }
  
  /**
   * Sets title property.
   * @param title                      title from the file.
   */
  public void setTitle(String title) { this.title = title; }
  
  /**
   * Sets the retention time of this Spectrum.
   * @param rt the retention time in seconds.
   */
  public void setRt(float rt) { this.rt = rt; }
 
  /**
   * Sets the fragmentation method of this Spectrum.
   * @param fragMethod the fragmentation method.
   */
  public void setActivationMethod(ActivationMethod fragMethod)
  {
	  this.activationMethod = fragMethod;
  }
  
  /**
   * Sets the ms level of this Spectrum.
   * @param msLevel the ms level.
   */
  public void setMsLevel(int msLevel)
  {
	  this.msLevel = msLevel;
  }
  
/****** FUNCTIONS *****/
  @Override
  public String toString() {
    return "Spectrum - mz: " +getPrecursorPeak().getMz() + ", peaks: " + size();
  }
  
  /**
   * Gets a clone of this spectrum with no peaks
   * @return a new spectrum without peaks
   */
  public Spectrum getCloneWithoutPeakList() {
    Spectrum newSpec = new Spectrum();
    newSpec.precursor = (Peak)this.precursor.clone();
    newSpec.startScanNum = this.startScanNum;
    newSpec.endScanNum = this.endScanNum;
    newSpec.title = this.title;
    newSpec.seqList = this.seqList;
    newSpec.annotation = this.annotation;
    newSpec.seqList = this.seqList;
    return newSpec;
  }
  
  /**
   * Gets a deep clone of this spectrum 
   * @return a new spectrum which is a copy of this
   */
  public Spectrum getDeepClone(){
	Spectrum newSpec = getCloneWithoutPeakList();
	for(Peak p : this){
		newSpec.add(p.clone());
	}
	return newSpec;
  }
  
  
  /**
   * Append a peak to the end of this spectrum.
   * @param peak                       peak object to be added.
   */
  public void addPeak(Peak peak) { 
    this.add(peak); 
  }


  /**
   * Correct parent mass according to annotation. Do nothing if annotation == null.
   */
  public void correctParentMass()
  {
	if(this.annotation == null || this.getCharge() <= 0)
		return;
	else
		this.precursor.setMz((annotation.getParentMass()+precursor.getCharge()*(float)Composition.PROTON)/precursor.getCharge());
  }
  
  /**
   * Correct parent mass according to parentMass.
   *  @param parentMass                
   */
  public void correctParentMass(float parentMass)
  {
	  this.precursor.setMz((parentMass+precursor.getCharge()*(float)Composition.PROTON)/precursor.getCharge());
  }
  
  
  /**
   * Correct parent mass according to input peptide. 
   */
  public void correctParentMass(Peptide pep)
  {
	if(this.getCharge() <= 0)
		return;
	else
		this.precursor.setMz((pep.getParentMass()+precursor.getCharge()*(float)Composition.PROTON)/precursor.getCharge());
  }
  
  /**
   * Correct charge of the precursor and set the charge of the peaks of this
   * spectrum to max(1, charge - 1).
   * @param charge                charge
   */
  public void setCharge(int charge){
	  this.precursor.setCharge(charge);
  }
  
  /**
   * Set the charge of the precursor peak
   * @param charge the new charge to use
   */
  public void setPrecursorCharge(int charge) {
    this.precursor.setCharge(charge);
  }
  
  /**
   * Returns a list of peaks that match the target mass within the tolerance
   * value. The absolute distance between the target mass and a returned peak
   * is less or equal that the tolerance value. The current implementation 
   * cycles through all peaks per call.
   * @param mass                       target mass.
   * @param tolerance               	tolerance.
   * @return an ArrayList object of the matching peaks. The array will be empty
   *         if there are no peaks within tolerance.
   */
  public ArrayList<Peak> getPeakListByMass(float mass, Tolerance tolerance) {
	  float toleranceDa = tolerance.getToleranceAsDa(mass);
	  return getPeakListByMassRange(mass-toleranceDa, mass+toleranceDa);
  }
  
  /**
   * Returns the most intense peak that is within tolerance of the target mass.
   * The current implementation takes linear time.
   * @param mass                       target mass.
   * @param tolerance                  tolerance.
   * @return a Peak object if there is match or null otherwise.
   */
  public Peak getPeakByMass(float mass, Tolerance tolerance)
  {
	  ArrayList<Peak> matchList = getPeakListByMass(mass, tolerance);
	  if(matchList == null || matchList.size() == 0)
		  return null;
	  else
		  return Collections.max(matchList, new IntensityComparator());
  }
  
  /**
   * Returns a list of peaks that match the target mass within the specified range.
   * Assuming spectrum is sorted by mass!!!
   * @param minMass                  minimum mass.
   * @param maxMass                  maximum mass.
   * @return an ArrayList object of the matching peaks. The array will be empty
   *         if there are no peaks within tolerance.
   */
  public ArrayList<Peak> getPeakListByMassRange(float minMass, float maxMass) {
	  ArrayList<Peak> matchList = new ArrayList<Peak>();
	    int start = Collections.binarySearch(this, new Peak(minMass,0,0));
	    if(start<0)
	    	start = -start-1;
	    for(int i = start; i<this.size(); i++)
	    {
	    	Peak p = this.get(i);
	    	if(p.getMz() > maxMass)
	    		break;
	    	else 
	    		matchList.add(p);
	    }
	  return matchList;
  }
  
  /**
   * Goes through the peaks and set their rank by intensity.
   * Rank 1: higest intensity peak
   */
  public void setRanksOfPeaks() {
    ArrayList<Peak> intensitySorted = new ArrayList<Peak>(this);
    Collections.sort(intensitySorted, Collections.reverseOrder(new IntensityComparator()));
    for(int i=0; i<intensitySorted.size(); i++) {
      intensitySorted.get(i).setRank(i+1);
    }
  }
  
  /**
   * Sets intensities of the charge two parent ion and its water loss to 0
   * @deprecated
   */
  public void filterPrecursorPeaks(Tolerance tolerance)
  {
	  filterPrecursorPeaks(tolerance, 0, 0);
  }
  
  /**
   * Filter (charge-reduced) precursor peaks with the specified offset
   */
  public void filterPrecursorPeaks(Tolerance tolerance, int reducedCharge, float offset)
  {
	  int c = this.getCharge() - reducedCharge;
	  float mass = (this.getParentMass()+c*(float)Composition.PROTON)/c+offset;
	  for(Peak p : getPeakListByMass(mass, tolerance))
		  p.setIntensity(0);
  }

  public void filterPrecursorPeaksAroundPM()
  {
	  for(int i=0; i<this.size(); i++)
	  {
		  float m = get(i).getMass();
		  int nominalMass = Math.round(m*Constants.INTEGER_MASS_SCALER);
		  if(nominalMass < 38)
			  this.get(i).setIntensity(0);
	  }
	  
	  // Remove all peaks with masses >= M+H - 38
	  int nominalPM = Math.round((getParentMass()-(float)Composition.H2O)*Constants.INTEGER_MASS_SCALER);
	  for(int i=this.size()-1; i>=0; i--)
	  {
		  float m = get(i).getMass();
		  int nominalMass = Math.round(m*Constants.INTEGER_MASS_SCALER);
		  if(nominalPM-nominalMass >= 38)
			  break;
		  this.get(i).setIntensity(0);
	  }

  }

  
  /**
   * The order of the a spectrum is determined by their parent masses.
   */
  public int compareTo(Spectrum s) {
    if(getParentMass() > s.getParentMass())
      return 1;
    else if(getParentMass() < s.getParentMass())
      return -1;
    return 0;
  }
  
  /*
  public void outputDtaFile(String fileName)
  {
    PrintStream out = null;
    try {
      out = new PrintStream(new File(fileName));
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    out.println(getParentMass()+Composition.H + " " + getCharge());
    for(Peak p : this)
    {
      out.println(p.getMass() + " " + p.getIntensity());
    }
    out.close();
  }
  
  public void outputPklFile(String fileName)
  {
    PrintStream out = null;
    try {
      out = new PrintStream(new File(fileName));
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    out.println(precursor.getMass() + " " + precursor.getIntensity() + " " + getCharge());
    for(Peak p : this)
    {
      out.println(p.getMass() + " " + p.getIntensity());
    }
    out.close();
  }

  public void outputPkl(PrintStream out)
  {
    out.println(precursor.getMass() + " " + precursor.getIntensity() + " " + getCharge());
    for(Peak p : this)
      out.println(p.getMass() + "\t" + p.getIntensity());
  }
  */

  /**
   * Output this spectrum to the input PrintStream as the mgf format.
   * It needs to be changed later.
   * @param out PrintStream object that the mgf spectrum will be written.
   */
  public void outputMgf(PrintStream out)
  {
    out.println("BEGIN IONS");
    if(this.title != null)
      out.println("TITLE=" + getTitle());
    if(this.annotation != null)
      out.println("SEQ=" + getAnnotationStr());
    
    float precursorMz = precursor.getMz();
    out.println("PEPMASS=" + precursorMz);
    if(startScanNum > 0)
      out.println("SCANS=" + startScanNum);
    int charge = getCharge();
    out.println("CHARGE=" + charge+ (charge>0 ? "+" : ""));
    for(Peak p : this)
    	if(p.getIntensity() > 0)
    		out.println(p.getMz() + "\t" + p.getIntensity());
    out.println("END IONS");
  }
  
  public void outputMgf(PrintStream out, boolean writeActivationMethod)
  {
          out.println("BEGIN IONS");
          if(this.title != null)
                  out.println("TITLE=" + getTitle());
          if(this.annotation != null)
                  out.println("SEQ=" + getAnnotationStr());
          if(this.getActivationMethod() != null && writeActivationMethod)
                  out.println("ACTIVATION=" + this.getActivationMethod().getName());
          float precursorMz = precursor.getMz();
          out.println("PEPMASS=" + precursorMz);
          if(startScanNum > 0)
                  out.println("SCANS=" + startScanNum);
          int charge = getCharge();
          out.println("CHARGE=" + charge+ (charge>0 ? "+" : ""));
          for(Peak p : this)
                  if(p.getIntensity() != 0)
                          out.println(p.getMz() + "\t" + p.getIntensity());
          out.println("END IONS");
  }

  
  
  // tmp function...
  public void outputMgf(PrintStream out, ArrayList<Peak> precursors)
  {
    out.println("BEGIN IONS");
    if(this.title != null)
      out.println("TITLE=" + getTitle());
    if(this.annotation != null)
      out.println("SEQ=" + getAnnotationStr());
    
   // float precursorMz = precursor.getMz();
	out.print("PEPMASS=");
    for(Peak p : precursors){
    	out.print(p.getMz()+",");
    }
    out.println();
    
    if(startScanNum > 0)
      out.println("SCANS=" + startScanNum);
    //int charge = getCharge();
    out.print("CHARGE=");//=" + charge+ (charge>0 ? "+" : "")
    for(Peak p : precursors){
    	out.print(p.getCharge()+ (p.getCharge()>0 ? "+" : "") +",");
    }
    out.println();

    for(Peak p : this)
    	if(p.getIntensity() > 0)
    		out.println(p.getMz() + "\t" + p.getIntensity());
    out.println("END IONS");
  }
  
  
  /**
   * Output this spectrum to the input PrintStream as the dta format.
   * It needs to be changed later.
   * @param fileName dta file name.
   */
  public void outputDta(String fileName)
  {
		PrintStream out = null; 
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(fileName)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
    out.println(this.getParentMass()+Composition.H+"\t"+this.getPrecursorPeak().getCharge());
    for(Peak p : this)
      out.println(p.getMz() + "\t" + p.getIntensity());
    out.close();
  }
  
  /**
   * Convert this spectrum into a dta string representation.
   * @return the dta representation.
   */
  public String toDta() {
    StringBuffer sb = new StringBuffer();
    sb.append(this.getParentMass()+Composition.H+"\t"+this.getPrecursorPeak().getCharge()+"\n");
    for(Peak p : this) sb.append(p.getMz() + "\t" + p.getIntensity() + "\n");
    return sb.toString();
  }
  
  /**
   * Inner class for intensity sorting.
   */	
  class IntensityComparator implements Comparator<Peak> {
    
    
    /**
     * Determines the order of peak objects, when sorted by intensity.
     * @param o1 the first peak object.
     * @param o2 the second peak object.
     * @return 1 if o1 > o2, -1 if o2 > o1 and 0 otherwise.
     */
    public int compare(Peak o1, Peak o2) {
      if(o1.getIntensity() > o2.getIntensity())  return 1;
      if(o2.getIntensity() > o1.getIntensity())  return -1;
      if(o1.getMz() > o2.getMz())            return 1;
      if(o2.getMz() > o1.getMz())            return -1;
      return 0;
    }
    
    
    /**
     * Determines equality according to intensities.
     * @param o1 the first peak object.
     * @param o2 the second peak object.
     * @return true if they compare to 0, false otherwise.
     */
    public boolean equals(Peak o1, Peak o2) {
      return compare(o1, o2) == 0;
    }
    
  }
  
  /*
   * Returns the most intense peak that is within tolerance of the target mass.
   * The current implementation takes linear time.
   * @deprecated
   * @param mass                       target mass.
   * @param toleranceDa                tolerance in Daltons.
   * @return a Peak object if there is match or null otherwise.
  public Peak getPeakByMass(float mass, float toleranceDa)
  {
    ArrayList<Peak> matchList = getPeakListByMassDa(mass, toleranceDa);
    if(matchList == null || matchList.size() == 0)
    	return null;
    else
    	return Collections.max(matchList, new IntensityComparator());
  }
	*/

  /*
   * @deprecated
   * @param mass
   * @param tolerancePPM
   * @return
  public Peak getPeakByMass(float mass, int tolerancePPM)
  {
    ArrayList<Peak> matchList = getPeakListByMassDa(mass, tolerancePPM);
    if(matchList == null || matchList.size() == 0)
    	return null;
    else
    	return Collections.max(matchList, new IntensityComparator());
  }
   */
  
}