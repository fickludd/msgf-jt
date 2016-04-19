package se.lth.immun.msgf.mzml;

import java.util.ArrayList;
import java.util.HashMap;

import se.lth.immun.msgf.MzMLCache;
import se.lth.immun.msgf.MzMLCache.SpecData;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;

public class MzMLSpectraMap implements SpectrumAccessorBySpecIndex {

	private MzMLCache mzmlCache = null;
	private HashMap<String, Integer> idMap = new HashMap<String, Integer>();
	
	public MzMLSpectraMap(MzMLCache mzmlCache) {
		this.mzmlCache = mzmlCache;
		for (int i = 0; i < mzmlCache.spectra().length(); i++) 
			idMap.put(mzmlCache.spectra().apply(i).id(), i);
	}

	public Spectrum getSpectrumBySpecIndex(int specIndex) {
		return ToSpectrum.from(specIndex, mzmlCache.spectra().apply(specIndex-1));
	}

	@Override
	public Spectrum getSpectrumById(String specId) {
		int specIndex = idMap.get(specId);
		return ToSpectrum.from(specIndex, mzmlCache.spectra().apply(specIndex-1));
	}
	
	public ArrayList<Integer> getSpecIndexList() {
		ArrayList<Integer> ret = new ArrayList<Integer>(mzmlCache.spectra().length());
		for (int i = 0; i < mzmlCache.spectra().length(); i++) ret.set(i, i);
		return ret;
	}

	@Override
	public String getID(int specIndex) {
		return mzmlCache.spectra().apply(specIndex-1).id();	
	}

	@Override
	public Float getPrecursorMz(int specIndex) {
		SpecData sd = mzmlCache.spectra().apply(specIndex-1);
		return (float)sd.precs()[0].precMz();
	}

	@Override
	public String getTitle(int specIndex) {
		return null;
	}
}