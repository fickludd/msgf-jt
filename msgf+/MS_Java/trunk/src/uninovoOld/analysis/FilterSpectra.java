package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;

public class FilterSpectra {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String s;
		BufferedLineReader in = new BufferedLineReader("/home/kwj/UniNovo/MSGFDBNoM16.txt");
		HashSet<Integer> sns = new HashSet<Integer>();
		while((s=in.readLine())!=null){
			if(s.startsWith("#")) continue;
			String[] token = s.split("\t");
			if(Float.parseFloat(token[14]) > 0.01f) continue;
			
			sns.add(Integer.parseInt(token[1]));
		}
		in.close();
		
		PrintStream out = new PrintStream("/home/kwj/workspace/inputs/CIDETDPairs/spectra/Trypsin/CID_filtered.mgf");
		Iterator<Spectrum> iterator = new SpectraIterator("/home/kwj/workspace/inputs/CIDETDPairs/spectra/Trypsin/CID.mgf", new MgfSpectrumParser());
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			
			if(sns.contains(spec.getScanNum())) continue;
			
			spec.outputMgf(out);
		}
		out.close();
	}

}
