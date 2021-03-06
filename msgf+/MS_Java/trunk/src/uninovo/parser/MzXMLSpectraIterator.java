package uninovo.parser;

import java.util.Iterator;

import uninovo.util.Spectrum;


/**
 * A data structure that allows iteration of the mzXML file.
 * Only MS/MS spectra will be returned.
 * @author jung
 *
 */
public class MzXMLSpectraIterator implements Iterator<Spectrum> {
	private Spectrum currentSpectrum;
	private boolean hasNext;
	private MzXMLSpectraMap map;
	private int scanNum;
	
	
	/**
	 * Constructor taking the file name.
	 * @param fileName
	 */
	public MzXMLSpectraIterator(String fileName) {
		this(fileName, 2, 2);
	}

	/**
	 * Constructor taking the file name and mslevel selectors
	 * @param fileName the path to the mzXML file.
	 * @param minMSLevel spectra with msLevel less than this will be ignored.
	 * @param maxMSLevel spectra with msLevel greater than this will be ignored.
	 */
	public MzXMLSpectraIterator(String fileName, int minMSLevel, int maxMSLevel) {
		map = new MzXMLSpectraMap(fileName).msLevel(minMSLevel, maxMSLevel);
		scanNum = 0;
		currentSpectrum = parseNextSpectrum();
		if(currentSpectrum != null)        hasNext = true;
		else                               hasNext = false;
	}
	
	/**
	 * Check whether there is more to parse.
	 * @return true if not done or false 
	 */
	@Override
	public boolean hasNext() {
		return hasNext;
	}
	
	
	/**
	 * Get next spectrum.
	 * @return the next spectrum.
	 */
	@Override
	public Spectrum next() {
		Spectrum curSpecCopy = currentSpectrum;
		currentSpectrum = parseNextSpectrum();
		if(currentSpectrum == null)
			hasNext = false;
		return curSpecCopy;
	}
	
	private Spectrum parseNextSpectrum()
	{
		Spectrum spec = null;
		while(++scanNum <= map.getScanCount())
		{
			spec = map.getSpectrumByScanNum(scanNum);
			if(spec != null)
				break;
		}
		return spec;
	}

	@Override
	public void remove() {
		assert(false);
	}	
}
