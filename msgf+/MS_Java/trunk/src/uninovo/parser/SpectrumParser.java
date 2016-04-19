package uninovo.parser;

import java.util.Hashtable;

import uninovo.util.Spectrum;

public interface SpectrumParser {
	Hashtable<Integer, Long> getScanNumMap(BufferedRandomAccessLineReader lineReader);
	public Spectrum readSpectrum(LineReader lineReader);
}
