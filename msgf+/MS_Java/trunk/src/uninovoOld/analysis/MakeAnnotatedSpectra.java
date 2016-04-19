package uninovoOld.analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;

import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;


public class MakeAnnotatedSpectra {
	public static void main(String[] args) throws IOException{
		String dir = "/home/kwj/workspace/inputs/CIDETDPairs/";
		String keyword = "Trypsin";
		String fragmentation = "CID";
		String fragmentation2 = "ETD";
		Boolean paired = true;
		float FDR = 0.01f;
		
		PrintStream out = new PrintStream(dir+"/annotated/"+fragmentation+ keyword + (paired? "_PAIRED" : "") + ".mgf");
		PrintStream out2 = null;
		if(paired) out2 = new PrintStream(dir+"/annotated/"+fragmentation2+ keyword + (paired? "_PAIRED" : "") + ".mgf");
		
		File[] list = new File(dir).listFiles();
		
		for(File file : list){
			if(!file.getAbsolutePath().contains(keyword)) continue;
			if(file.getAbsolutePath().endsWith((paired? "SUM.txt" : fragmentation+".txt"))){
				System.out.println(file.getAbsolutePath());
				BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
				String s, filename = "";
				HashMap<Integer, Float> scanNumPrecursorMz = new HashMap<Integer, Float>();
				HashMap<Integer, String> annotations = new HashMap<Integer, String>();
				
				while((s = in.readLine())!=null){
					if(s.startsWith("#")||s.isEmpty()) continue;
					String[] t = s.split("\t");
					filename = t[0];
					if(paired && filename.contains("novpredictAD")) continue;
					if(Float.parseFloat(t[t.length-1]) > FDR) continue;
					
					String[] index = t[1].split("/");
					//String[] frag = t[3].split("/");
					if(paired && index.length ==1){
						System.out.println("shit*************");
						continue;
					}
					for(int i=0; i<index.length; i++){
						//if(!fragmentation.equals(frag[i])) continue;
						scanNumPrecursorMz.put(Integer.parseInt(index[i]), Float.parseFloat(t[4]));
						annotations.put(Integer.parseInt(index[i]), t[7].substring(t[7].indexOf('.')+1, t[7].lastIndexOf('.')));
					}
				}
				in.close();
				System.out.println(filename + "\t" + scanNumPrecursorMz.size() + "\n" + dir + "spectra/" + filename);
			
				if(!filename.isEmpty()){
					Iterator<Spectrum> iterator = null;
					int scanoffset = 0;
					if(filename.endsWith(".mgf")){
						scanoffset = 1;
						iterator = new SpectraIterator(dir + "spectra/" + filename, new MgfSpectrumParser());
					}else{
						iterator = new MzXMLSpectraIterator(dir + "spectra/" + filename);
					}
					
					while(iterator.hasNext()){
						Spectrum spec = iterator.next();
						if(!scanNumPrecursorMz.containsKey(spec.getScanNum()+scanoffset))continue;
						
						if(spec.getActivationMethod().toString().equals(fragmentation)) {
							assert(scanNumPrecursorMz.get(spec.getScanNum()+scanoffset) == spec.getPrecursorPeak().getMz());
							spec.setAnnotation(new Peptide(annotations.get(spec.getScanNum()+scanoffset)));
							spec.outputMgf(out, true);
						}
						if(paired && spec.getActivationMethod().toString().equals(fragmentation2)){
							assert(scanNumPrecursorMz.get(spec.getScanNum()+scanoffset) == spec.getPrecursorPeak().getMz());
							spec.setAnnotation(new Peptide(annotations.get(spec.getScanNum()+scanoffset)));
							spec.outputMgf(out2, true);
						}
						
						
					}
				}
			}
		}
		out.close(); out2.close();
	}
}
