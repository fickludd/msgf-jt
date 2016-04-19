package kyowon_local;

import java.io.IOException;

import parser.BufferedLineReader;

public class countGPLen {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		//HashMap<Integer, Integer> gpnum = new HashMap<Integer, Integer>();
		//HashMap<Integer, Float> gplen = new HashMap<Integer, Float>();
		
		for(int len = 8;len <= 20; len++){
			int num = 0, lennum = 0;
			BufferedLineReader in = new BufferedLineReader("/home/kwj/workspace/outputs/shewLengthAll.grc");
			
			String s;
			
			boolean count = false;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#")){
					String[] token = s.split("\t");
					if(token[6].length() == len) count = true;
					else count = false;
				}else if(count){
					num++;
					String[] tmp = s.split(",");
					lennum += tmp.length;
				}
				
			}
			
			System.out.println(len + "\t:"+(float)lennum/num);
			in.close();
		}

	}

}
