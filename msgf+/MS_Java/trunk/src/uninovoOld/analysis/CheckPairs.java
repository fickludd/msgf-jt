package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MgfSpectrumParser;

public class CheckPairs {

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		String CID = "/home/kwj/workspace/inputs/Zubarev/spectra/mergedHCD2.mgf", 
				ETD = "/home/kwj/workspace/inputs/Zubarev/spectra/mergedETD.mgf";
		Iterator<Spectrum> i1 = new SpectraIterator(CID,
				new MgfSpectrumParser());
		Iterator<Spectrum> i2 = new SpectraIterator(ETD,
				new MgfSpectrumParser());

		PrintStream out1 = new PrintStream(CID + ".checked");
		PrintStream out2 = new PrintStream(ETD + ".checked");

		String prevFile1 = null, prevFile2 = null;
		boolean go1 = true, go2 = true;
		Spectrum s1 = null, s2 = null;
		int sn = 0;
		
		while (i1.hasNext() && i2.hasNext()) {
			
			if(go1) s1 = i1.next();
			if(go2) s2 = i2.next();

			String title1 = s1.getTitle(); // 1
			String title2 = s2.getTitle(); // 2

			String[] t1 = title1.split(":");
			String[] t2 = title2.split(":");

			if (prevFile1 == null) {
				prevFile1 = t1[0];
			}

			if (prevFile2 == null) {
				prevFile2 = t2[0];
			}

			String file1 = t1[0];
			String file2 = t2[0];

			if (!file1.equals(file2)) {
				if (file1.equals(prevFile1)) { // file1 is old
					go2 = false;
					go1 = true;
				} else if (file2.equals(prevFile2)) {
					go1 = false;
					go2 = true;
				}
			}else{
				go1 = go2 = true;
			}

			prevFile1 = file1;
			prevFile2 = file2;

			if(!go1 || !go2) continue;
			
			int diff = Integer.parseInt(t2[t2.length-1].trim()) - Integer.parseInt(t1[t1.length-1].trim());
			
			if(diff > 1){
				go2 = false;
				go1 = true;
				System.out.println(s1.getTitle() + "\t" +s2.getTitle());
			}else if(diff < 1){
				go1 = false;
				go2 = true;
				System.out.println(s1.getTitle() + "\t" +s2.getTitle());
			}else{
				go1 = go2 = true;
			}
	
			if(!go1 || !go2) continue;
			assert(s1.getParentMass() == s2.getParentMass());
			s1.outputMgf(out1);
			s2.outputMgf(out2);
			sn++;
			if(sn > 20000) break;
			//System.out.println(sn);
		}

		out1.close();
		out2.close();

	}

}
