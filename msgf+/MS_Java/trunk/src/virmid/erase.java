package virmid;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class erase {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/Users/kwj/Documents/workspace/virmid/joint");
		String outdir = "/Users/kwj/Documents/workspace/virmid/joint/fixed/";
		int i = 0;
		
		/*PrintStream outsh = new PrintStream(outdir + "st.sh");
		outsh.println("source ~/.bashrc");
		for(int j=1;j<=50;j++){
			outsh.println("qsub -l h_vmem=5G -v PATH=/data/home/kyjeong/python/bin/:/data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/:$PATH -v PYTHONPATH=/data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/joint_snv_mix/:$PYTHONPATH st"+j+".sh");
			outsh.println("sleep 5");
		}
		outsh.close();
		*/
		for(File file : dir.listFiles()){
			//if(file.getName().endsWith("st.sh")){
			//	continue;
			//}
			if(file.getName().endsWith(".sh")){
				BufferedLineReader in = new BufferedLineReader(file.getAbsolutePath());
				PrintStream out = new PrintStream(outdir + file.getName());
				String s;
				while((s=in.readLine())!=null){
					if(s.contains("-01A-") || s.contains("-10A-") || s.contains("c52847b9d5bb9da3b067e9d68f3b3a01") || s.contains("6bf16f633d53c28c9746de992fb66c0b")){
						System.out.println(++i + s);
						
						String[] t = s.split(" ");
						String newS = "";
						for(int z=0;z<=2;z++)
							newS += t[z]+" ";
						newS += t[4]+ " " + t[3] + " ";
						for(int z=5;z<t.length;z++)
							newS += t[z] + " ";
						s = newS;
					}
					//s = s.replace("--input_file:normal", "TTTTT");
					//s = s.replace("--input_file:tumor", "--input_file:normal");
					//s = s.replace("TTTTT", "--input_file:tumor");
				/*	s = s.replace("-01A-", "-TTTTTTT-");
					s = s.replace("-10A-", "-01A-");
					s = s.replace("-TTTTTTT-", "-10A-");
					s = s.replace("c52847b9d5bb9da3b067e9d68f3b3a01", "-XXXXX-");
					s = s.replace("6bf16f633d53c28c9746de992fb66c0b", "c52847b9d5bb9da3b067e9d68f3b3a01");
					s = s.replace("-XXXXX-", "6bf16f633d53c28c9746de992fb66c0b");
					*/
					out.println(s);
				}
				out.close();
				in.close();
			}
		}
	}

}
