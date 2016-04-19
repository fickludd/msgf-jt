package kyowon_local;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class MakeTrainingSpec {
	static public void main(String[] args) throws IOException{
		String inputmgf = "/home/kwj/workspace/inputs/Training/Zubarev_HCD_Annotated.mgf";
		PrintStream out = new PrintStream(inputmgf.substring(0, inputmgf.lastIndexOf(".")) + "_train.mgf");
		
		HashSet<String> peps = new HashSet<String>();
		Iterator<Spectrum> iterator = new SpectraIterator(inputmgf, new MgfSpectrumParser());
		while(iterator.hasNext()){
			Spectrum s= iterator.next();
			if(peps.contains(s.getAnnotationStr())) continue;
			peps.add(s.getAnnotationStr());
			s.outputMgf(out);
			
		}
		out.close();
	}
}
