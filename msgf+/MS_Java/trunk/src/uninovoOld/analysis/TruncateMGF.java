package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.Random;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class TruncateMGF {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String specFile = "/home/kwj/Dropbox/Training/CIDTrypsin_train.mgf";
		PrintStream out = new PrintStream(specFile+"truncated.mgf");
		
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		
		int i=0;
		Random r = new Random();
		
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			if(s.getCharge()!=2) continue;
			
			if(r.nextFloat()>1f/3) continue;
			
			s.outputMgf(out);
			
			
			
			if(i++>5000) break;
			
			
			
		}
		System.out.println(i);
		out.close();
	}

}
