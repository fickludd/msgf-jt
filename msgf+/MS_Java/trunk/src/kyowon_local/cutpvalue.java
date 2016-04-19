package kyowon_local;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import parser.BufferedLineReader;

public class cutpvalue {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		float t = 1e-9f;
		BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/ostrich/NEWOSTRICHMERGED.txt");
		BufferedWriter out = new BufferedWriter(new FileWriter("/home/kwj/workspace/outputs/ostrich/NEWOSTRICHMERGED.txt.after"));
	
		int filecol = 0;
		int speccol = 1;
		int pepcol = 5;
		int proteincol = 7;
		//int scorecol = 7;
		int pvalcol = 6;
		
		String[] trexpeps = {
			//	"GLVGAPGLRGLPGK",// no
				"GVVGLPGQR",
				"GVQGPPGPQGPR",//
				"GATGAPGIAGAPGFPGAR",//
				"GLPGESGAVGPAGPIGSR", //
				"GSAGPPGATGFPGAAGR",//
				"GAPGPQGPSGAPGPK",//
				"VNVADCGAEALAR"//
		};
		
		String[] bernpeps = {
			"GLAGPQGPR",//
			"GPPGESGAAGPTGPIGSR", 
			"GEPGPAGLPGPAGER",//
			"GPPGSSGSTGK", 
			"GAAGLPGVAGAPGLPGPR",//
			"PGCPGPMGEK", 
			"LSDLHAQK",//
	
			"RNVADCGAEALAR"
		};
		
		String[] sangtaepeps = {
				".NVADCGAEALAR",//
				"DEVTPA.VVVAR",
				"LVNELTEFAK",//
				"SSN.LSGSTLR", // 
				".GVDAGAAGDPER",//
				"S.IHVALVTGGNK",
				".LVNELTEFAK",
				"L.NELTEFAK",
				"ENAGEDPGLAR",
				"EDCLSG.KPK",//
				".GC.GPMGEK", 
			};
		
		
		HashSet<String> filescannums = new HashSet<String>();
		HashSet<String> peps = new HashSet<String>();
		HashMap<String, String> u = new HashMap<String, String>();;
		HashMap<String, HashSet<String>> species = new HashMap<String, HashSet<String>>();
		String s;
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#")) continue;
			String[] token = s.split("\t");
			
		//	if(!s.contains("Tyrannosaurus")){
				if(Float.parseFloat(token[pvalcol]) > t) continue;
				if(token[proteincol].contains("Trypsin")) continue;
				if(token[proteincol].contains("Keratin")) continue;
		//	}
			if(filescannums.contains(token[filecol]+token[speccol])) continue;
			filescannums.add(token[filecol]+token[speccol]);
			String pep = token[pepcol].substring(token[pepcol].indexOf('.')+1, token[pepcol].lastIndexOf('.'));
			peps.add(pep);
			
			String color = "black";
			
		//	for(String tp : sangtaepeps){
		//		if(pep.toUpperCase().matches(tp)) color = "orange";
				//else if(tp.matches(pep.toUpperCase())) color = "cyan";
		//	}
			
		
			for(String tp : bernpeps){
				if(pep.equalsIgnoreCase(tp)) color = "blue";
				else if(tp.contains(pep.toUpperCase())) color = "yellow";
			}
			
			for(String tp : trexpeps){
				if(pep.equalsIgnoreCase(tp)) color = "red";
				else if(tp.contains(pep.toUpperCase())) color = "green";
			}
		//	if(color.equals("black")) continue;

			out.write(s+"\n");
			
			if(color.equals("black")) continue;

			
			String key =token[proteincol].substring(token[proteincol].indexOf(' ')+1, token[proteincol].indexOf(" ", token[proteincol].indexOf("OS=")))+
			"&{\\color{"+color+"}{"+token[pepcol]+"}}&"+"&"+String.format("%.2e", Float.parseFloat(token[pvalcol]));

			if(!u.containsKey(pep)) u.put(pep, key);
			else{
				if(!(Float.parseFloat(u.get(pep).split("&")[3]) < Float.parseFloat(token[pvalcol])))
					u.put(pep, key);
			}

			if(!species.containsKey(pep)) species.put(pep, new HashSet<String>());
			HashSet<String> sp = species.get(pep);
			//System.out.println(token[proteincol]);
			sp.add(token[proteincol].substring(token[proteincol].indexOf("OS="), token[proteincol].indexOf(" ", token[proteincol].indexOf("OS="))));
			
		}
			
		//{\color{blue}{http://www.devdaily.com/}}
		in.close();out.close();
		System.out.println(filescannums.size() + " " + peps.size()+" " + u.size());
		
		float colnum =0, hemnum = 0, othernum=0;
		
		for(String pep : u.keySet()){
			String x = u.get(pep);
			System.out.println(x+"\\\\\\hline");//
			
			if(x.contains("Collagen")) colnum++;
			else if(x.contains("Hemoglobin")) hemnum++;
			else othernum++;
		}
		
		System.out.println(colnum+ " " + hemnum + " " + othernum);
		
		for(String pep : species.keySet()){
			System.out.println(pep + "\t" + species.get(pep));
		}
		
	}
	
}
