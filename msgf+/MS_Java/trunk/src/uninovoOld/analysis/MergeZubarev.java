package uninovoOld.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import msutil.ActivationMethod;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class MergeZubarev {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String ZubarevDir = "/home/kwj/workspace/inputs/Zubarev/";
		File dir = new File(ZubarevDir+"original");
		PrintStream out = null;
		
		for(File file : dir.listFiles()){
			if(!file.getAbsolutePath().endsWith("mgf")) continue;
			if(file.getAbsolutePath().contains("out")) continue;

			out = new PrintStream(ZubarevDir + file.getName());
		
			
			Iterator<Spectrum> iterator = new SpectraIterator(file.getAbsolutePath(), new MgfSpectrumParser());
			while(iterator.hasNext()){
				Spectrum s = iterator.next();
				if(s.getTitle().contains("experiment: 2")){
					s.setActivationMethod(ActivationMethod.HCD);//HCD
				}else{
					s.setActivationMethod(ActivationMethod.ETD);
				}
				s.outputMgf(out, true);
				
			}
			out.close();
		}
		
	}

}
