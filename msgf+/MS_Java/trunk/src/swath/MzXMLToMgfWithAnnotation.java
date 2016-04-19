package swath;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import msutil.Peptide;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MzXMLSpectraIterator;

public class MzXMLToMgfWithAnnotation {

	
	
	private static HashMap<String, HashMap<Integer,String>> getAnnotations(String msgfFileName){
		HashMap<String, HashMap<Integer,String>> annotations = new HashMap<String, HashMap<Integer,String>>();
		File msgfFile = new File(msgfFileName);
		File[] msgfFileList;
		if(msgfFile.isDirectory())
			msgfFileList = msgfFile.listFiles(new misc.FileFilter.FileExtFilter("txt"));
		else
		{
			msgfFileList = new File[1];
			msgfFileList[0] = msgfFile;
		}
		
		for(File file : msgfFileList){
			String s;
			try {
				BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
				while((s=in.readLine())!=null){
					if(s.startsWith("#")) continue;
					String[] token = s.split("\t");
					
					String fn = token[0];
					if(Float.parseFloat(token[13]) > 0.01) continue;
					int sn = Integer.parseInt(token[1]);
					String pep = token[7].substring(token[7].indexOf('.')+1, token[7].lastIndexOf('.'));
					
					if(!annotations.containsKey(fn)) annotations.put(fn, new HashMap<Integer, String>());
					HashMap<Integer, String> v = annotations.get(fn);
					v.put(sn, pep + ":" + token[6]);
				}
				in.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			
			
		}
		return annotations;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		
		String msgf = "/home/kwj/workspace/inputs/SWATH/new";
		String mzxml = "/home/kwj/workspace/inputs/SWATH/new";
		HashMap<String, HashMap<Integer,String>> annotations = getAnnotations(msgf);
		try {
			PrintStream out;
			
			out = new PrintStream(mzxml+ System.getProperty("file.separator") + "converted.mgf");
			
			for(String file : annotations.keySet()){
				String fn = mzxml + System.getProperty("file.separator") + file;
				if(new File(fn).exists()){
						
						HashMap<Integer,String> k = annotations.get(file);
						System.out.println(fn);
						MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fn);
						while(itr.hasNext())
						{
							Spectrum spec = itr.next();
							String pep = k.get(spec.getScanNum());
							if(pep != null){
								String[] t = pep.split(":");
								spec.setAnnotation(new Peptide(t[0]));
								spec.setCharge(Integer.parseInt(t[1]));
								spec.outputMgf(out);
							}
						}
						
					
				}
			}
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
