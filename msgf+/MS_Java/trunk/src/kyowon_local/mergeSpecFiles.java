package kyowon_local;

import java.io.IOException;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class mergeSpecFiles {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Iterator<Spectrum> iterator = null;
		String spectrumFileName = null;
		
		try {
			iterator = new SpectraIterator(spectrumFileName, new MgfSpectrumParser());
			
			
		}catch  (IOException e){
			
		}

	}

}
