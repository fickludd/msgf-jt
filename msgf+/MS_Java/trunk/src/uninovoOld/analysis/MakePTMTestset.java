package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;

import msutil.Enzyme;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;

public class MakePTMTestset {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String specFile="/home/kwj/Dropbox/Test/Modified/modified_spectra_nontrivial_nonambiguous_.05_FLR.mgf";
		String annotationFile = "/home/kwj/Dropbox/Test/Modified/annotations.txt";
		PrintStream out = new PrintStream("/home/kwj/Dropbox/Test/Modified/test.mgf");
		
		int charge = 2;
		HashSet<String> annotations = new HashSet<String>();
		
		BufferedLineReader in = new BufferedLineReader(annotationFile);
		String s;
		while((s=in.readLine())!=null){
			annotations.add(s);
		}
		in.close();
		int sn = 0;
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(spec.getCharge() != charge || !annotations.contains(spec.getTitle())) continue;
			
			String a = spec.getTitle();
			System.out.println(a);
			
		
			a=a.replace("(", "");
			a=a.replace(")", "");
			a=a.replace(",", "");
			if(!a.contains("-18.0106")){
				a=a.replace("S-18", "S+79.966");
				a=a.replace("T-18", "T+79.966");
				a=a.replace("Y-18", "Y+79.966");
			}else{
				a=a.replace("S-18.0106", "S+79.966");
				a=a.replace("T-18.0106", "T+79.966");
				a=a.replace("Y-18.0106", "Y+79.966");
			}
			
			
			System.out.println(a);
			
			Peptide annotation = new Peptide(a);
			
			System.out.println(annotation + "\t" + annotation.getParentMass() + "\t" + spec.getParentMass());
			if(Math.abs(annotation.getParentMass() - spec.getParentMass()) > 2.0) continue;
			spec.setAnnotation(annotation);
			if(!Enzyme.TRYPSIN.isCleaved(annotation)) continue;
			spec.setTitle(a);
			//spec.correctParentMass();
			//
			spec.outputMgf(out);
			sn++;
		}
		
		in = new BufferedLineReader("/home/kwj/Dropbox/Test/Modified/test.mgf");
		out = new PrintStream("/home/kwj/Dropbox/Test/Modified/test.mgf2");
		String ann = "";
		while((s=in.readLine())!=null){
			if(s.startsWith("TITLE")) ann = s.split("=")[1];
			else if(s.startsWith("SEQ")) {
				String t = "SEQ=" + ann;
				s = t;
			}
			out.println(s);
		}
		out.close();
		in.close();
		 
		System.out.println(sn);
		out.close();
	}

}
