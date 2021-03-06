package uninovo.parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import uninovo.util.Pair;
import uninovo.util.Spectrum;
import uninovo.util.SpectrumAccessorByScanNum;


/**
 * A data structure that allows iteration of the spectrum file by increasing order of parent mass.
 * @author sangtae
 *
 */
public class SortedSpectraIterator implements Iterator<Spectrum> {
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();
		String fileName = "/home/sangtaekim/Research/Data/HeckWhole/Spectra/090121_NM_Trypsin_20.mzXML";
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fileName);
		MzXMLSpectraMap map = new MzXMLSpectraMap(fileName);
		SortedSpectraIterator sortedItr = new SortedSpectraIterator(itr, map);
		int numSpecs = 0;
		while(sortedItr.hasNext())
		{
			Spectrum spec = sortedItr.next();
			System.out.println(spec.getParentMass());
			numSpecs++;
		}
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
	}
	private Spectrum currentSpectrum;
	private boolean hasNext;
	private int index = -1;
	private SpectrumAccessorByScanNum map;
	private final int numSpecs;
	
	private ArrayList<Integer> scanNumList;

	
	/**
	 * Constructor taking the file name.
	 * @param fileName
	 */
	public SortedSpectraIterator(Iterator<Spectrum> itr, SpectrumAccessorByScanNum map) {
		int numSpecs = 0;
		ArrayList<Pair<Integer,Float>> scanNumPMPairList = new ArrayList<Pair<Integer,Float>>();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			scanNumPMPairList.add(new Pair<Integer,Float>(spec.getScanNum(), spec.getParentMass()));
			numSpecs++;
		}
		this.numSpecs = numSpecs;
		Collections.sort(scanNumPMPairList, new Pair.PairComparator<Integer,Float>(true));
		scanNumList = new ArrayList<Integer>();
		for(Pair<Integer,Float> p : scanNumPMPairList)
			scanNumList.add(p.getFirst());
		
		this.map = map;
		index = -1;
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
		++index;
		if(index >= scanNumList.size())
			return null;
		int scanNum = scanNumList.get(index);
		return map.getSpectrumByScanNum(scanNum);
	}

	@Override
	public void remove() {
		assert(false);
	}	
	
	/**
	 * Returns the number of spectra.
	 * @return the number of spectra.
	 */
	public int size() {
		return numSpecs;
	}
}
