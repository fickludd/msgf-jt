package swath;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import parser.BufferedLineReader;

public class EvalProjectedCosSecond {

	/**
	 * @param args
	 * @throws IOException 
	 */
	
	
	
	static HashMap<Integer, Float> getCosTable(String file, String delm, String proteinKey) throws IOException{
		BufferedLineReader in = new BufferedLineReader(file);
		String s;
		HashMap<Integer, Float> out = new HashMap<Integer, Float>();
		
		while((s=in.readLine())!=null){
			if(s.startsWith("Parsed") || s.startsWith("#") || s.startsWith("matching") || s.contains("DECOY")) continue;
			s = s.replace("Scan Number: ", "ScanNumber:");
			String[] token = s.split(delm);
			if(token.length < 8)
				System.out.println(s);
			
			if(!proteinKey.isEmpty() && !token[token.length-1].contains(proteinKey)) continue;
			if(token[token.length-1].toLowerCase().contains("keratin")) continue;
			
			int libsn = Integer.parseInt(token[8].split(":")[1].trim());
			if(libsn == Integer.parseInt(token[1])) continue;
			
			float cos = Float.parseFloat(token[7]);
			if(out.containsKey(libsn)){
				float prevcos = out.get(libsn);
				if(prevcos < cos) out.put(libsn, cos);
			}else 
				out.put(libsn, cos);
		}
		
		in.close();
		
		return out;
	}
	
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		String ida = "/home/kwj/workspace/inputs/SWATH/Ecoli/ACG_swathdevelopment_14341_UPS1-400fm_UPS1_vs_IDA_library_projcos_min0.6.txt.proteinAdded.txt";
		String swath = "/home/kwj/workspace/inputs/SWATH/Ecoli/ACG_swathdevelopment_14342_UPS1-400fm_SWATH_5600_vs_IDA_library_projcos_min0.6.txt.proteinAdded.txt";
		
		String proteinKey = "coli";
		
		String out = ida+ "_" + proteinKey + ".m";
		
		PrintStream outs = new PrintStream(out);
		HashMap<Integer, Float> scosval = getCosTable(swath, "\t", proteinKey);
		HashMap<Integer, Float> icosval = getCosTable(ida, "\t", proteinKey);
		
		outs.println("a = [");
		for(int k : scosval.keySet()){
			if(icosval.containsKey(k)){
				outs.println(scosval.get(k) + "\t" + icosval.get(k) + "\t" + k);
			}
			
		}
		outs.println("];");
		outs.close();
		
	}

}
