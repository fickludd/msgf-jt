package kyowon_local;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import parser.BufferedLineReader;

public class mm {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		
		File folder = new File("/home/kwj/workspace/inputs/HEK/");
		File[] files = folder.listFiles();
		BufferedWriter out = new BufferedWriter(new FileWriter("/home/kwj/workspace/inputs/HEK/result.mgf"));
		for(File file : files){
			int n=0;
			BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
			String s;
			
			while((s=in.readLine()) != null && n <= 150){
				out.write(s+"\n");
				if(s.startsWith("END IONS")) n++;
			}
			in.close();
		}
		
		out.close();
	}

}
