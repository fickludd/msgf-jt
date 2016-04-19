package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class SeparateZubarev {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String input = "/home/kwj/Dropbox/Test/Blind/19Feb_ADfulldigest_HCDETDTop5_2h_f0_01.mgf";
		String output = "/home/kwj/Dropbox/Test/Blind/19Feb_ADfulldigest_HCDETDTop5_2h_f0_01_ETD.mgf";
		String keyWord = "experiment: 1";
		
		PrintStream out = new PrintStream(output);
		Iterator<Spectrum> iterator = new SpectraIterator(input, new MgfSpectrumParser());

		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			
			if(!s.getTitle().contains(keyWord)) continue;
			s.outputMgf(out);
		}
		out.close();
	}

}
