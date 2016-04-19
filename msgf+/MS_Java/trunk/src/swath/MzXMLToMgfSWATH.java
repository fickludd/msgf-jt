package swath;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import msgf.Tolerance;
import msutil.Peak;
import msutil.Spectrum;
import parser.MzXMLSpectraIterator;
import parser.PSM;
import parser.PSMList;
import parser.SWATHParser;

public class MzXMLToMgfSWATH {
	private static HashMap<Integer,PSM> getPSMs(String swathFileName){
		HashSet<Integer> snsWithMultipleAnnotations = new HashSet<Integer>();
		HashMap<Integer,PSM> ret = new HashMap<Integer, PSM>();
		
		PSMList<PSM> psms = SWATHParser.parse(swathFileName);
		for(PSM psm : psms){
			int sn = psm.getScanNum();
			
			if(ret.containsKey(sn)){
				if(!ret.get(sn).getPeptideStr().equals(psm.getPeptideStr())){
					snsWithMultipleAnnotations.add(sn);
				}
			}
			ret.put(sn, psm);
		}
	//	System.out.println(ret.size());
		for(int sn : snsWithMultipleAnnotations)
			ret.remove(sn);
	//	System.out.println(ret.size());
		psms = new PSMList<PSM>(ret.values());
		
		psms = psms.getDistinctiveSpectralSet();//.getDistinctivePeptideSet();
		ret.clear();
		
		for(PSM psm : psms){
			int sn = psm.getScanNum();
			ret.put(sn, psm);
		}
		
		return ret;
	}
	
	public static void main(String[] args) throws IOException {
		String fn = "/home/kwj/workspace/inputs/SWATH/sw/cdk47ct_swath-cdk47ct.mzXML";
		String libfn = "/home/kwj/workspace/inputs/SWATH/pw/cdk47ct_IDA-cdk47ct.mzXML";
				
		String txt = "/home/kwj/workspace/inputs/SWATH/sw/gringar_SLRIACG01_cdk47ct_swath_msplit_libmsgfdbids_IDs.txt";
		PrintStream out = new PrintStream("/home/kwj/workspace/inputs/SWATH/sw/converted.mgf");
		PrintStream outlib = new PrintStream("/home/kwj/workspace/inputs/SWATH/sw/convertedlib.mgf");
		
		Tolerance tol = new Tolerance(96.08f, true);
		
		HashMap<Integer,PSM> psms = getPSMs(txt);
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fn);
		HashMap<String, Spectrum> os = new HashMap<String, Spectrum>();
		HashMap<String, Float> ic = new HashMap<String, Float>();
		HashMap<String, Integer> libsn = new HashMap<String, Integer>();
		HashMap<Integer, Spectrum> libss = new HashMap<Integer, Spectrum>();
		
		while(itr.hasNext()){
			Spectrum s = itr.next();
			if(psms.containsKey(s.getScanNum())){
				PSM psm = psms.get(s.getScanNum());
				s.setAnnotation(psm.getPeptide());
				s.setCharge(psm.getCharge());
				
				float ionc = 0;
				for(Peak p : s){
					ionc += p.getIntensity();
				}
				
				if(os.containsKey(s.getAnnotationStr())){
					
					if(ic.get(s.getAnnotationStr()) < ionc){
						ic.put(s.getAnnotationStr(), ionc);
						os.put(s.getAnnotationStr(), s);
						libsn.put(s.getAnnotationStr(), (int)psm.getScore("LSN"));
					}
				}else{
					ic.put(s.getAnnotationStr(), ionc);
					os.put(s.getAnnotationStr(), s);
					libsn.put(s.getAnnotationStr(), (int)psm.getScore("LSN"));
				}
			}
		}
		
		//15345
		//12183
		
		if(!libfn.isEmpty()){
			itr = new MzXMLSpectraIterator(libfn);
			while(itr.hasNext()){
				Spectrum s = itr.next();
				if(libsn.containsValue(s.getScanNum())){
					libss.put(s.getScanNum(), s);
				}
			}
		}
		
		System.out.println(os.size());
		for(String t : os.keySet()){
			Spectrum s = os.get(t);
			Spectrum so = s.getCloneWithoutPeakList();
			for(Peak p : s){
				Peak po = p.clone();
				po.setMz(p.getMz() + tol.getToleranceAsDa(p.getMz()));
				so.add(po);
			}
			so.outputMgf(out);
			if(!libfn.isEmpty()){
				Spectrum ls = libss.get(libsn.get(t));
				if(ls != null) ls.outputMgf(outlib);
			}
			
			
		}
		outlib.close();
		out.close();
	}
}
