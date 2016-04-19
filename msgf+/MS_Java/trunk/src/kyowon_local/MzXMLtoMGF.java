package kyowon_local;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import msutil.Spectrum;
import parser.MzXMLSpectraIterator;

public class MzXMLtoMGF {
	static public void main(String[] args) throws IOException{
		Iterator<Spectrum> iterator = null;
		
		iterator = new MzXMLSpectraIterator(args[0]);
		PrintStream out = new PrintStream(args[1]);
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			spec.outputMgf(out);
		}
		
		out.close();
	}
}
