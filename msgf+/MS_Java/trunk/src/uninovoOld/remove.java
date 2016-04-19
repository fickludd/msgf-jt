package uninovoOld;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class remove {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String in = "/home/kwj/workspace/inputs/Zubarev/annotated/HCD_PAIRED.mgf";
		String in2 = "/home/kwj/workspace/inputs/Zubarev/annotated/ETD_PAIRED.mgf";
		
		String ref ="/home/kwj/workspace/inputs/Training/HCD_train.mgf";
		String out = "/home/kwj/Dropbox/Test/HCDETD/HCD.mgf";
		String out2 = "/home/kwj/Dropbox/Test/HCDETD/ETD.mgf";
		
		Iterator<Spectrum> iterator = new SpectraIterator(ref, new MgfSpectrumParser());
		HashSet<String> annos = new HashSet<String>();
		
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			
			annos.add(s.getAnnotationStr());
		}
		
		iterator = new SpectraIterator(in, new MgfSpectrumParser());
		Iterator<Spectrum> iterator2 = new SpectraIterator(in2, new MgfSpectrumParser());
		
		PrintStream outmgf = new PrintStream(out);
		PrintStream outmgf2 = new PrintStream(out2);
		int[] num = new int[4];
		
		while(iterator.hasNext()){
			Spectrum s = iterator.next();
			Spectrum s2 = iterator2.next();
			if(s.getCharge() > 3 || s.getCharge() < 2) continue;
			
			if(s.getAnnotation().isModified() || s.getAnnotation().size() < 5) continue;
			if(!s.getAnnotation().equals(s2.getAnnotation())){
				System.out.print("ss");
				System.exit(0);
			}
			if(annos.contains(s.getAnnotationStr())) continue;
			//if(num[s.getCharge()] ++ >= 1000) continue;
			
			annos.add(s.getAnnotationStr());
			s.outputMgf(outmgf);
			s2.outputMgf(outmgf2);
			
		}
		
		outmgf.close();
		outmgf2.close();
	}

}
