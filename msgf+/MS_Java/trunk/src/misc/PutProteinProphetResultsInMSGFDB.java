package misc;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import parser.BufferedLineReader;

public class PutProteinProphetResultsInMSGFDB {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String msgfResult = args[0];//"/home/kwj/workspace/ProteoSafeNew/worker/tasks/0d16ddc911724fee8f7ce7e2e7ef3e05/finalResult/cdaa4b76acc241f99d597a46e427766a";//args[0];
		String proteinProphetResult = args[1];//"/home/kwj/workspace/ProteoSafeNew/worker/tasks/0d16ddc911724fee8f7ce7e2e7ef3e05/msgfProtXMLResult/bdaedfc288bd4422a8c095b77229a7a4.pep.inter.xml";//args[1];
		String msgfResultProteinProphetResultAppended = args[2];//"/home/kwj/workspace/ProteoSafeNew/worker/tasks/0d16ddc911724fee8f7ce7e2e7ef3e05/finalResult/cdaa4b76acc241f99d597a46e427766a2";//args[2];
		
		HashMap<String, String> map = new HashMap<String, String>();
		BufferedLineReader in = new BufferedLineReader(proteinProphetResult);
		String s;
		boolean set = false;
		String pep = null;
		String prot = null;
		while((s=in.readLine())!=null){
			if(s.startsWith("<spectrum_query")){
				set = true;						
			}			
			if(!set) continue;
			
			if(s.startsWith("<search_hit")){
				int t = s.indexOf("protein=\"") + 9;		
				char pre = s.charAt(s.indexOf("peptide_prev_aa=\"") + 17);
				char post = s.charAt(s.indexOf("peptide_next_aa=\"") + 17);
				prot = s.substring(t, s.indexOf('\"', t)) + "(pre=" + pre + ",post=" + post + ")";	
				
				t = s.indexOf("peptide=\"") + 9;
				pep = s.substring(t, s.indexOf('\"', t));
			}
			
			if(pep != null && prot != null){
				map.put(pep, prot);
			}
				
			if(s.startsWith("</spectrum_query>")){
				set = false;
				prot = null;
				pep = null;
			}
		}
		
		in.close();
		
		BufferedLineReader in2 = new BufferedLineReader(msgfResult);
		PrintStream out = new PrintStream(msgfResultProteinProphetResultAppended);
		
		while((s=in2.readLine())!=null){
			if(s.startsWith("#")){
				out.println(s+"\tProteinProphetProtein");
				continue;
			}
			String[] token = s.split("\t");
			pep = token[8];
			boolean mapped = false;
			if(map.containsKey(pep)){
				String protein = map.get(pep);
				if(protein != null){
					mapped = true;
					s = s + "\t" + protein;
				}else{
					s = s + "\tN/A";
				}
			}else s = s + "\tN/A";
			if(!mapped) System.out.println("Matching PSM not found : " + pep);
			
			out.println(s);
		}
		
		in2.close();
		out.close();
		
	}

}
