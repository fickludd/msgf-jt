package hybridsearch;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import parser.BufferedLineReader;

public class MergeTargetDecoy {

	
	
	private static void updateBestScoringPSMs(String file, int snCol, int scoreCol, boolean isGreaterBetter, HashMap<String, String> bestScoringPSMs){
		BufferedLineReader in;
		try {
			in = new BufferedLineReader(file);
			String s;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#")|| s.startsWith("Scan#")){
					bestScoringPSMs.put("HEADER", s);
					continue;
				}
				
				String[] token = s.split("\t");
				
				String key = token[snCol];
				
				float score = Float.parseFloat(token[scoreCol]);
				boolean isCurrentBetter = true;
				
				if(bestScoringPSMs.containsKey(key)){
					float prevScore = Float.parseFloat(bestScoringPSMs.get(key).split("\t")[scoreCol]);
					isCurrentBetter = isGreaterBetter ? score > prevScore : score < prevScore;			
				}
				
				if(isCurrentBetter) bestScoringPSMs.put(key, s);
				
			}
			
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	
	
	public static void main(String[] args) throws IOException {
		String target = args[0];
		String decoy = args[1];
		boolean isGreaterBetter = args[2].equals("1");
		int snCol = Integer.parseInt(args[3]);
		int scoreCol = Integer.parseInt(args[4]);
		
		
		PrintStream out = new PrintStream(args[5]);
		
		HashMap<String, String> bestScoringPSMs = new HashMap<String, String>();
		
		updateBestScoringPSMs(target, snCol, scoreCol, isGreaterBetter, bestScoringPSMs);
		updateBestScoringPSMs(decoy, snCol, scoreCol, isGreaterBetter, bestScoringPSMs);
		
		if(bestScoringPSMs.containsKey("HEADER"))
			out.println(bestScoringPSMs.get("HEADER"));
		
		for(String key : bestScoringPSMs.keySet()){
			if(key.equals("HEADER")) continue;
			
			out.println(bestScoringPSMs.get(key));
		}
		
		out.close();
	}

}
