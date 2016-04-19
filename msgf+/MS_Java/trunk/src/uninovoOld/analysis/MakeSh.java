package uninovoOld.analysis;

import java.io.File;

public class MakeSh {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String Dir = "/home/kwj/workspace/inputs/Zubarev";
		File dir = new File(Dir);
		int mode = 4;
		int enzyme = 1;

		for(File file : dir.listFiles()){
			if(!file.toString().endsWith("mgf")) continue;
			String t = "CID";
			if(mode == 2) t = "ETD";
			if(mode == 3) t= "HCD";
			if(mode == 4) t="SUM";
			String s =
				"java -jar -Xmx2000M MSGFDB.jar -s " + file.getAbsolutePath() + 
				" -d /home/kwj/workspace/inputs/DataBases/ipi.HUMAN.v3.73.fasta -t 20ppm -o " + 
				file.getAbsolutePath() + t + ".txt" + " -tda 1 -c13 0 -m "+ mode +" -e "+enzyme+" -mod /home/kwj/workspace/MSGFDB/Mods.txt -maxCharge 5";
			System.out.println(s);
		}
	}

}
