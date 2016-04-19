package uninovoOld.analysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import parser.BufferedLineReader;

public class AnalyzeTrex {
	public static void main(String[] args) throws IOException{
		String file = "/home/kwj/UniNovo/ETDUni.txtmod";
		float FDR = .01f;
		boolean isSpectrumFDR = true;
		
		HashMap<String, String> tpep = new HashMap<String, String>();
		HashSet<String> dpep = new HashSet<String>();
		
		HashSet<Integer> tspec = new HashSet<Integer>();
		HashSet<Integer> dspec = new HashSet<Integer>();
		 
		float fdr = 100f;//0.008757157	2969	1.9753085E-13	spec num : 4101	qspec num : 4442

		float th = 1e-8f;//3.8415167E-14f;
		
		int sn = 0;
		int tsn = 0;
		int[] histo = null;
		while(fdr > FDR){			
			BufferedLineReader in = new BufferedLineReader(file);
			String s;
			
			tpep.clear();
			dpep.clear();
			tspec.clear();
			dspec.clear();
			
			int pn = 0;
			sn = 0;
			histo = new int[4000];
			String[] k = null;
			while((s=in.readLine())!=null){
				if(s.startsWith("#")) continue;
				String[] token = s.split("\t");

				int cn = Integer.parseInt(token[1]);
				if(Float.parseFloat(token[6])>th) continue;
				
				if(pn!=cn){
					sn++;
					if(k!=null){
						//if(Float.parseFloat(k[11])>=1)
						{
							String pep = k[5].substring(k[5].indexOf('.'), k[5].lastIndexOf('.'));
							if(k[8].startsWith("REV_")){
								dpep.add(pep);
								dspec.add(Integer.parseInt(k[1]));
							}
							else{ 
								//
								tpep.put(pep, k[3]+"\t"+ k[11]);
								tspec.add(Integer.parseInt(k[1]));
								histo[Math.round(Float.parseFloat(k[11]))+200]++;
							}
						}
					}
					k = token;
					
					
				}else if(k!=null){
					if(Float.parseFloat(token[7]) > Float.parseFloat(k[7])){
						k = token;
					}else if(Float.parseFloat(token[7]) == Float.parseFloat(k[7]) && Float.parseFloat(token[11]) < Float.parseFloat(k[11])){
						k = token;
					}
				}
				
				
				pn = cn;
				
				
			}
			
			if(tsn == 0) tsn = sn;
			
			fdr = (float)dpep.size()/tpep.size();
			if(isSpectrumFDR) fdr = (float)dspec.size()/tspec.size();
			
			System.out.println(fdr+"\t"+tspec.size());
			if(fdr > 2f* FDR)
				th /= 2f;
			else th /= 1.5f;
			
		}
		
		for(String pep : tpep.keySet())
			System.out.println(pep + "\t" + tpep.get(pep));
		
		for(int i=0;i<histo.length;i++){
			if(histo[i]>60) System.out.println(i-200 + "\t" + histo[i]);
		}
		
		System.out.println((float)dpep.size()/tpep.size() + "\t" + tpep.size() + "\t" + th + "\tspec num : " + sn + "\tqspec num : " + tsn);
		
		
		
		/*
		String MSGFDBoutfile = "/home/kwj/workspace/inputs/Heck_DDDT/CIDH.out";
		String s;
		HashSet<String> identifiedByMSGFDB = new HashSet<String>();
		
		file= "/home/kwj/workspace/inputs/Heck_DDDT/CID_HighRes_Tryp.mgf";
		PrintStream out = new PrintStream("/home/kwj/workspace/inputs/Heck_DDDT/CID_HighRes_Tryp_filtered.mgf");
		
		BufferedLineReader in = new BufferedLineReader(MSGFDBoutfile);
		while((s=in.readLine())!=null){
			if(s.startsWith("#")) continue;
			String[] token = s.split("\t");
			
			if(Float.parseFloat(token[14]) > FDR) continue;
			
			identifiedByMSGFDB.add(token[1]+"\t"+token[4]);//TODO what happens for CID/ETD?? check...
		}
		in.close();
		
		Iterator<Spectrum> iterator = new SpectraIterator(file, new MgfSpectrumParser());
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(identifiedByMSGFDB.contains((spec.getScanNum()+1)+"\t"+spec.getPrecursorPeak().getMz()))continue;
				
			spec.outputMgf(out);
		}
		*/
	}
}
