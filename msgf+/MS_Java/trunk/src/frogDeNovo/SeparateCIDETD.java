package frogDeNovo;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class SeparateCIDETD {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String f = "/home/kwj/workspace/inputs/Nebedev/Nobel prize_Lessonae-ox-zip1-no-excl.mgf";
		
		PrintStream C = new PrintStream(f+"CID.mgf");
		PrintStream E = new PrintStream(f+"ETD.mgf");

		Iterator<Spectrum> iH = new SpectraIterator(f, new MgfSpectrumParser());
		
		while(iH.hasNext()){
			Spectrum s = iH.next();
			if(s.getTitle().contains("CID")){
				s.outputMgf(C);
			}else s.outputMgf(E);
		}
		
		C.close();
		E.close();
	}

}
