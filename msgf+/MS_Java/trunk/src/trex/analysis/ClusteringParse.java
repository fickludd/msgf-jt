package trex.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import parser.BufferedLineReader;

public class ClusteringParse {
	static public void main(String[] args) throws IOException{
		ArrayList<String> clusters = new ArrayList<String>();
		
		File folder = new File("/home/kwj/workspace/inputs/ostrich/clust");
		File[] files = folder.listFiles();
		
		for(File file : files){
			if(!file.getAbsolutePath().endsWith(".clust")) continue;
			
			BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
			
			String s;
			String c="";
			
			while((s=in.readLine())!=null){
				if(s.isEmpty()){
					String[] t = c.split("\n");
					
					boolean towrite = false;
					boolean towrite2 = false;
					for(int i=1; i<t.length; i++){
						
						if(Integer.parseInt(t[i].split("\t")[1]) == 15){
							towrite = true;
						}else{
							towrite2 = true;
						}
					}
					
					if(towrite && towrite2) clusters.add(c);
					
					c="";
				}else				
					c+=s+"\n";
				
			}
			
			in.close();
			
		}
		
		for(String cluster : clusters){
			System.out.println("******************");
			System.out.println(cluster);
		}
	}
}
