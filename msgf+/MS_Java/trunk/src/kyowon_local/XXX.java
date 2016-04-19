package kyowon_local;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import parser.BufferedLineReader;

public class XXX {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File folder = new File(args[0]);
		
		for(File file : folder.listFiles()){
			if(!file.getName().endsWith(".txt")) continue;
			///data/kyjeong/outputs/inspectSearchResults/Decoy
			BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
			BufferedWriter out = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".rev"));
			
			String s;
			
			while((s=in.readLine())!=null){
				if(s.startsWith("#")){
					out.write(s+"\n");
					continue;
				}
				
				String[] token = s.split("\t");
				
				for(int i=0;i<token.length-1;i++){
					if(i!=3)
						out.write(token[i]+"\t");
					else 
						out.write("XXX_"+token[i]+"\t");
				}
				
				out.write(token[token.length-1]+"\n");
				
			}
			
			out.close();
			in.close();
		}
		
		
		
		
	}

}
