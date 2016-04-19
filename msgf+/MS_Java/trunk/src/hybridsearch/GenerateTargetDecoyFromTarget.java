package hybridsearch;

import java.io.IOException;

import msdbsearch.ReverseDB;

public class GenerateTargetDecoyFromTarget {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String in = null;
		String out = null;
		boolean genDecoy = true;
		
		for(int i=0; i<args.length; i++){
			String arg = args[i];
			
			if(arg.equals("-i")) in = args[i+1];
			if(arg.equals("-o")) out = args[i+1];
			if(arg.equals("-t")) genDecoy = args[i+1].equals("1");
		}
		
		//System.out.println(in);

		//System.out.println(out);

		//System.out.println(key);
		try {		
			if(genDecoy){
				ReverseDB.reverseDB(in, in+"_rev.fasta");
			}
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}

}
