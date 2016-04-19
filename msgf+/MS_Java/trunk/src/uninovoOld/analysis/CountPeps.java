package uninovoOld.analysis;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Peptide;
import parser.BufferedLineReader;

public class CountPeps {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String name = "HCD.txt";//1564 vs 1508 hcd , 3422 vs 3400 cid , 2906 vs 2809 etd ,  3549 vs 3546 cidetd
		String name1 = "iPRG.txt";
		String file = "/home/kwj/Dropbox/MSGF/"+name;// 850 vs 825 hcdetd
		String file2 = "/home/kwj/Dropbox/" + name1 + "mut";//score mut 287 / 1 : mod 412 / 4 
		String file3 = "/home/kwj/Dropbox/" + name1 + "mod2";//score mut 287 / 1 : mod 412 / 4 
		String grc = "/home/kwj/Dropbox/" + name1.substring(0, name1.indexOf("txt")) + "grc";
		
		String out = "/home/kwj/Dropbox/" + name1 + "out";
		// 06	342	3	0.877193 when filter out
		Enzyme enzyme = null;//Enzyme.TRYPSIN;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		BufferedLineReader in = new BufferedLineReader(file);
		BufferedLineReader in2 = new BufferedLineReader(file2);
		BufferedLineReader in3 = new BufferedLineReader(file3);
		BufferedLineReader gin = new BufferedLineReader(grc);
		
		//ETYGEMADC+126CAK
		//System.out.println( new Peptide("ETYGEMADCCAK").getMass());//717.77
		
		int maxCharge = 10;
		boolean specNonOverlap = false;
		boolean pepNonOverlap = false;
		
		int pepCol = 7; // 2
		int specProbCol = 11;
		int chargeCol = 6;
		int snCol=2;
		float specProb = 1.151E-9f;// 1.151E-9f; HCD
		
		int pepCol2 = 5; // 2
		int specProbCol2 = 12;
		int scoreCol = 7;
		int chargeCol2 = 4;
		int snCol2=1;

		float score = -20f;//60.6f;// 130.1f
		float specProb2 = 2e-10f;

		int[] PTMs = new int[401];
		
		HashSet<Integer> denovosns = new HashSet<Integer>();
		HashSet<Integer> sns = new HashSet<Integer>();
		HashSet<Integer> sns2 = new HashSet<Integer>();
		HashSet<Integer> sns2t = new HashSet<Integer>();
		HashSet<Integer> sns2d = new HashSet<Integer>();
		HashMap<Integer, String> sns2StringMap = new HashMap<Integer, String>();
		HashSet<String> peps = new HashSet<String>();
		HashSet<String> pepsd = new HashSet<String>();
		HashSet<String> peps2 = new HashSet<String>();
		HashSet<String> peps2d = new HashSet<String>();
		HashMap<String, String> pepAnnotationMap = new HashMap<String, String>();
		HashMap<String, Integer> mutationCountMap = new HashMap<String, Integer>();
		
		String s;
		
		while((s=gin.readLine())!=null){
			if(s.startsWith("#")){
				int sn = Integer.parseInt(s.split("\t")[snCol2]);
				denovosns.add(sn);
			}
		}
		
		/*
		while((s=in.readLine())!=null){
			if(s.startsWith("#")) continue;
			String[] token = s.split("\t");
			if(Integer.parseInt(token[chargeCol]) > maxCharge) continue;
			if(Float.parseFloat(token[specProbCol]) > specProb) continue;
			String anno = token[pepCol];
			String pep = anno.substring(anno.indexOf('.')+1, anno.lastIndexOf('.')); 
			
			if(enzyme!=null && !(anno.charAt(anno.length()-1) == '_' || enzyme.isCleaved(new Peptide(pep)) && (anno.charAt(0) == '_' || enzyme.isCleavable(anno.charAt(0))))){
				continue;
			}
			
			for(String n : token[snCol].split("/"))
				sns.add(Integer.parseInt(n));
			
			if(!s.contains("REV_")){
				peps.add(pep);
			}else{
				if(!peps.contains(pep))
					pepsd.add(pep);
			}
		}
		
		while((s=in2.readLine())!=null){
			if(s.startsWith("#")) continue;
			String[] token = s.split("\t");
			if(Integer.parseInt(token[chargeCol2]) > maxCharge) continue;	
			int sn = Integer.parseInt(token[snCol2]);
			if(specNonOverlap && sns.contains(sn)) continue;
			
			String anno = token[pepCol2];
			String pep = anno.substring(anno.indexOf('.')+1, anno.lastIndexOf('.'));
			if(pepNonOverlap && (peps.contains(pep) || pepsd.contains(pep))) continue;
			
			sns2StringMap.put(sn, s);
			
		}*/
		
		while((s=in3.readLine())!=null){
			if(s.startsWith("#")) continue;
			if(s.equals("null")) continue;
			String[] token = s.split("\t");
			//System.out.println(s);
			if(Integer.parseInt(token[chargeCol2]) > maxCharge || Integer.parseInt(token[chargeCol2]) == 0) continue;	
			int sn = Integer.parseInt(token[snCol2]);
			if(specNonOverlap && sns.contains(sn)) continue;
			
			String anno = token[pepCol2];
			String pep = anno.substring(anno.indexOf('.')+1, anno.lastIndexOf('.'));
			if(pepNonOverlap && (peps.contains(pep) || pepsd.contains(pep))) continue;
			
			if(sns2StringMap.containsKey(sn)){
				//float ss = Float.parseFloat(sns2StringMap.get(sn).split("\t")[scoreCol]);
				//float thisss = Float.parseFloat(token[scoreCol]);
				float ss = Float.parseFloat(sns2StringMap.get(sn).split("\t")[sns2StringMap.get(sn).split("\t").length-1]);
				float thisss = Float.parseFloat(token[token.length-1]);
				
				/*Peptide p = new Peptide(pep);
				for(AminoAcid aa : p){
					
					int ptm = Math.round(aa.getMass() - aaSet.getAminoAcid(Character.toUpperCase(aa.getResidue())).getMass());
					if(ptm == 22){
						System.out.println(sns2.get(sn));
					}
				}*/
				
				if(thisss < ss) sns2StringMap.put(sn, s);
			}else
				sns2StringMap.put(sn, s);
			
		}
		
		HashMap<Integer, String> output = new HashMap<Integer, String>();
		
		for(int sn : sns2StringMap.keySet()){
			s = sns2StringMap.get(sn);
			String[] token = s.split("\t");
	
			if(Float.parseFloat(token[token.length-1]) > specProb2) continue;
			if(Float.parseFloat(token[scoreCol]) < score) continue;

			String anno = token[pepCol2];
			String pep = anno.substring(anno.indexOf('.')+1, anno.lastIndexOf('.'));
			
			if(enzyme!=null && !(anno.charAt(anno.length()-1) == '*' || enzyme.isCleaved(new Peptide(pep)) && (anno.charAt(0) == '*' || enzyme.isCleavable(anno.charAt(0))))){
				continue;
			}
			
			sns2.add(sn);
			
			//if(peps2.contains(pep)) continue;
			
			if(!s.contains("REV_")){
				sns2t.add(sn);
				peps2.add(pep);
				output.put(sn, s);
				if(pepAnnotationMap.containsKey(pep)){

					float ps = Float.parseFloat(pepAnnotationMap.get(pep).split("\t")[pepAnnotationMap.get(pep).split("\t").length-1]);
					if(ps > Float.parseFloat(s.split("\t")[s.split("\t").length-1]))

						pepAnnotationMap.put(pep, s);
				}else pepAnnotationMap.put(pep, s);
				
				
				
				if(token.length < 14){
					Peptide p = new Peptide(pep);
					for(AminoAcid aa : p){
						
						int ptm = Math.round(aa.getMass() - aaSet.getAminoAcid(Character.toUpperCase(aa.getResidue())).getMass());
						if(ptm!=0){
						//	if(aa.getResidue() == 'c')// C modification
						//		ptm += 57;//
							PTMs[200 + ptm]++;
							break;
						}
					}
				}else{
					String mut = token[15].substring(token[15].indexOf("[")+1, token[15].indexOf("]"));
					int c = 0;
					
					String[] tt = mut.split("->");
					AminoAcid aa1 = aaSet.getAminoAcid(tt[0].charAt(0));
					AminoAcid aa2 = aaSet.getAminoAcid(tt[1].charAt(0));
					
					float mass1 = aa1 == null? 0 : aa1.getMass();
					float mass2 = aa2 == null? 0 : aa2.getMass();
					
					mut += "\t" + (mass2 - mass1);
					if(mutationCountMap.containsKey(mut)){
						c = mutationCountMap.get(mut);
					}
					mutationCountMap.put(mut, c+1);
				}
			}else{
				sns2d.add(sn);
				if(!peps2.contains(pep))
					peps2d.add(pep);
			}
		}
		
		PrintStream outfile = new PrintStream(out);
		
		int maxsn = 0;
		for(int i : output.keySet()){
			maxsn = Math.max(maxsn, i);
		}
		
		for(int i=0; i<=maxsn; i++){
			if(output.containsKey(i)){
				outfile.println(output.get(i));
			}
		}
		
		outfile.close();
		int c = 0;
		for(String pep : peps2){
			if(peps.contains(pep)){
				c++;
			}else{
				String[] t = pepAnnotationMap.get(pep).split("\t");
				if(pep.contains("+") || pep.contains("-"))
					if(Math.round(Float.parseFloat(t[t.length-1]))==-48)
					System.out.println("PTM: " + pep + "\t" + t[chargeCol2] + "\t" + t[snCol2] + "\t" + t[pepCol2] + "\t" + t[t.length-1]);
			}
		}
		
		for(String pep : peps2){
			if(peps.contains(pep)){
			}else{
				String[] t = pepAnnotationMap.get(pep).split("\t");
				if(!pep.contains("+") && !pep.contains("-"))
					System.out.println("MUT: " + pep + "\t" + t[chargeCol2] + "\t" + t[snCol2] + "\t" + t[pepCol2] + "\t" + t[t.length-1]);
			}
		}
		
		
		int d=0;
		for(int sn : sns){
			if(sns2StringMap.containsKey(sn)){
				d++;
			}else{
		//		System.out.println(pep);
			}
		}
		
		
		//System.out.println(c+"\n" + sns.size() + "\t" + sns2.size() + "\t" + d+ "\n" +peps2d);
		

		int modnum = 0;
		System.out.println("PTM=[");
		for(int i = 0; i<PTMs.length; i++){
			System.out.println((i-200) + "\t" + PTMs[i]);
			modnum += PTMs[i];
		}
		System.out.println("];");

	
		int mutnum = 0;
		for(String mut : mutationCountMap.keySet()){
			System.out.println(mut + "\t" + mutationCountMap.get(mut));
			mutnum += mutationCountMap.get(mut);
		}
		
		HashSet<Integer> tsns = new HashSet<Integer>();		
		tsns.addAll(sns);
		tsns.addAll(sns2);
		HashSet<String> tpeps = new HashSet<String>();
		tpeps.addAll(peps);
		tpeps.addAll(peps2);
		HashSet<String> tpepsd = new HashSet<String>();
		tpepsd.addAll(pepsd);
		tpepsd.addAll(peps2d);
		
		System.out.println("Mod Num: " + modnum + "\tMut Num: " + mutnum);
		
		System.out.println("\nTool\tSpec#\tPep#\tDecoyPep#\tFDR");
		System.out.println("MSGFDB: " + sns.size() + "\t" +peps.size()+"\t" + pepsd.size() + "\t" + ((float)pepsd.size() /peps.size() *100));
		System.out.println("AdaNovo: " + sns2.size() + "\t"+ peps2.size()+"\t" + peps2d.size() + "\t" + ((float)peps2d.size() /peps2.size() *100)+"\t" + ((float)sns2d.size() /sns2t.size() *100));
		System.out.println("Total: " + tsns.size() + "\t" + (tpeps.size())+"\t" + (tpepsd.size()) + "\t" + ((float)(tpepsd.size())/(tpeps.size())  *100));
	
		int remainingsn = denovosns.size();
		for(int sn : denovosns){
			if(tsns.contains(sn)){
				remainingsn--;
			}
		}
		System.out.println("Remaining high quality spectra: " + remainingsn);
		//
		in.close();
		in2.close();
	}

}
