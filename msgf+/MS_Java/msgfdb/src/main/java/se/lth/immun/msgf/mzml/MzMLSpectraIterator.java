package se.lth.immun.msgf.mzml;

import java.util.Iterator;

import edu.ucsd.msjava.msutil.Spectrum;

import se.lth.immun.msgf.MzMLCache;

public class MzMLSpectraIterator implements 
		Iterator<edu.ucsd.msjava.msutil.Spectrum>, 
		Iterable<edu.ucsd.msjava.msutil.Spectrum> 
{

	private MzMLCache mzmlCache = null;

	private int i = 0;
	
	public MzMLSpectraIterator(MzMLCache mzmlCache) {
		this.mzmlCache = mzmlCache;
	}
	
	public boolean hasNext() {
		return i < mzmlCache.spectra().length();
	}

	/**
	 * Get next spectrum.
	 * @return the next spectrum.
	 */
	public Spectrum next() {
		int index = i+1;
		return ToSpectrum.from(index, mzmlCache.spectra().apply(i++));
	}
	
	public Iterator<Spectrum> iterator() {
		return this;
	}

	public void remove() {
		throw new UnsupportedOperationException("MzMLSpectraIterator.remove() not implemented");
	}
}