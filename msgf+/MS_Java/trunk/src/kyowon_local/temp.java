package kyowon_local;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import msutil.AminoAcidSet;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MS2SpectrumParser;

public class temp {
	static public void main(String[] argv) throws IOException{
		BufferedLineReader fasta = new BufferedLineReader("/home/kwj/workspace/inputs/Sequest/ids.fasta");
		SpectraIterator iterator = new SpectraIterator("/home/kwj/workspace/inputs/Sequest/selected.ms2", new MS2SpectrumParser());
		PrintStream out = new PrintStream("/home/kwj/workspace/inputs/Sequest/selected.mgf");
		ArrayList<Peptide> ids = new ArrayList<Peptide>();
		
		String s;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		while((s=fasta.readLine())!=null){
			if(s.startsWith(">")) continue;
			ids.add(new Peptide(s.substring(1, s.length()-1), aaSet));
		}
		
		fasta.close();
		int i=0;
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			spec.setAnnotation(ids.get(i++));
			spec.correctParentMass();
			spec.outputMgf(out);
		}
		
		out.flush();out.close();
	}
}
