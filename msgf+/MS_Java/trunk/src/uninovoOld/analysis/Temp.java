package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class Temp {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String training = "/home/kwj/Dropbox/Training/HCD_train.mgf";
		String spec1 = "/home/kwj/workspace/inputs/Zubarev/annotated/HCD_PAIRED.mgf";
		String spec2 = "/home/kwj/workspace/inputs/Zubarev/annotated/ETD_PAIRED.mgf";
		int specNum = 1000;
		int charge = 3;
		
		PrintStream out1 = new PrintStream(spec1 + "Paired" + charge + ".mgf");
		PrintStream out2 = new PrintStream(spec2 + "Paired" + charge + ".mgf");
		
		HashSet<String> titlesInTraining = new HashSet<String>();
		
		Iterator<Spectrum> iterator = new SpectraIterator(training, new MgfSpectrumParser());
		
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			titlesInTraining.add(s.getTitle());
		}
		
		iterator = new SpectraIterator(spec1, new MgfSpectrumParser());
		Iterator<Spectrum> iterator2 = new SpectraIterator(spec2, new MgfSpectrumParser());
		
		HashSet<String> peps = new HashSet<String>();
		int sn = 0;
		
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			Spectrum t = iterator2.next();
			
			if(s.getCharge() != charge || s.getAnnotation().isModified() || s.getAnnotation().size()<5) continue;
			if(titlesInTraining.contains(s.getTitle())) continue;
			if(peps.contains(s.getAnnotationStr())) continue;
			
			
			peps.add(s.getAnnotationStr());
			s.outputMgf(out1);
			t.outputMgf(out2);
			
			sn ++;
			
			if(sn >= specNum) break;
		}
		System.out.println(sn);
		out1.close();
		out2.close();

	}

}
