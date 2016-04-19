package virmid;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import parser.BufferedLineReader;

public class TestSpecificity {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		File dir = new File("/home/kwj/workspace/MutectSpecificity");
		
		for(File file : dir.listFiles()){
			if(file.isDirectory()) continue;
			String fn = file.getAbsolutePath();
			//if(!fn.endsWith("post.txt")) continue;

			BufferedLineReader in = new BufferedLineReader(fn);
			String s;
			
			int cnt=0;
			while((s=in.readLine())!=null){
				if(s.startsWith("chrom"))continue;
				if(s.startsWith("##")||s.startsWith("contig"))continue;
				
				if(!s.endsWith("REJECT")){
					cnt++;
				}
				/*String[] t = s.split("\t");
				//String ch = tok[0];
				float p=1-Float.parseFloat(t[8])-Float.parseFloat(t[12])-Float.parseFloat(t[16]);
				if(p>=0.5f){
					cnt++;
				}*/
			}
			System.out.println(file.getName()+": " + cnt);
			
			in.close();
			
		}
		
	}

}
