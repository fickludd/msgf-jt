package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class tmp {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String specFile = "/home/kwj/workspace/inputs/Zubarev/spectra/mergedETDchecked.mgf";
		PrintStream out = new PrintStream("/home/kwj/workspace/inputs/Zubarev/spectra/mergedETDchecked2.mgf");
		
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		
		int i=0;
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			if(i == 2){
				i=0;
				s.outputMgf(out);
			}
			i++;
			
			
		}
		out.close();
	}

}
