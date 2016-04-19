package tdaAnalysis;

import java.io.File;

public class ShGenerator {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String folder = "/home/kwj/workspace/inputs/TDA/Kyowon/";
		
		for(int i=3;i<=3;i++){
			
			String tol = "2.5Da -c13 0";
			if(i == 28 || (i >=16 && i<=18)) 
				tol = "30ppm -c13 1";
			
			int k=i;
			if(i==28) k = 27;
			if(!new File(folder + k).exists()) continue;
			
			File mgfFolder = new File(folder + k + "/mgf/");
			File dbFolder = new File(folder + k + "/db/");
			
			for(File mgf : mgfFolder.listFiles()){
				int j=0;
				if(!mgf.getName().endsWith("mgf")) continue;
				for(File db : dbFolder.listFiles()){
					if(db.getName().endsWith("fasta")){
						System.out.println("java -Xmx2000m -jar MSGFDB.jar -s " + mgf + " -d " + db + " -t " + tol + " -o " + folder + k + "/" + j + (i==28 ? "_high" : "") + ".txt"   
							+" -tda 0 -e 1");
						j++;
					}
				}
			}
			
		}
	}

}
