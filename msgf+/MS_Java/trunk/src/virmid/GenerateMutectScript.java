package virmid;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class GenerateMutectScript {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String header = "#!/bin/bash\n#\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n#";
		String backbone1 = "java -Xmx2g -jar /data/home/kyjeong/Mutect/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence /data/home/kyjeong/reference/Homo_sapiens_assembly19.fasta --input_file:normal ";
		String backbone2 = " --input_file:tumor ";
		//String backbone3 = "---out";
		
		String data = "/home/kwj/Dropbox/Virmid/sensitivity/ourrecal.txt";
		String trainsh = "/home/kwj/workspace/sh/nmt";
		String exeshFinal = "/home/kwj/workspace/sh/nmt.sh";
		
		//String runsh = "/home/kwj/workspace/sh/rt";
		//String exeshFinal2 = "/home/kwj/workspace/sh/rt.sh";
	
		
		PrintStream exeshFinalo = new PrintStream(exeshFinal);
		//PrintStream exeshFinalo2 = new PrintStream(exeshFinal2);
		exeshFinalo.println("source ~/.bashrc");
		BufferedLineReader datain = new BufferedLineReader(data);
		String s;
		//make -C /data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/results/result1
		//sho.println(header);
		int cntr = 1;
		File dir = new File("/home/kwj/workspace/MutectCCMS/sh/");
		
		while((s=datain.readLine())!=null){
			String[] t = s.split(" ");
			
			boolean finished = false;
			for(File file : dir.listFiles()){
				if(file.getName().contains("mt"+cntr+".sh.o")){
					BufferedLineReader tmp = new BufferedLineReader(file.getAbsolutePath());
					String sss;
					while((sss=tmp.readLine())!=null){
						if(sss.contains("NSRuntimeProfile")){
							finished = true;
							break;
						}
					}
					tmp.close();
				}
			}
			if(finished){
				cntr++;
				continue;
			}
			
			PrintStream exesho = new PrintStream(trainsh + cntr + ".sh");
			//PrintStream exesho2= new PrintStream(runsh + cntr + ".sh");
			
			exesho.println(header);
			exesho.print(backbone1+t[0]+backbone2 + t[1] + " --out /data/home/kyjeong/Mutect/results2" + t[0].substring(t[0].lastIndexOf("/")) + ".txt");
	
			exesho.println(" --coverage_file /data/home/kyjeong/Mutect/results2" + t[0].substring(t[0].lastIndexOf("/")) + ".cov.txt");
			
			exeshFinalo.println("qsub -l h_vmem=5G nmt"+cntr + ".sh\nsleep 10");
			
			
			
			//exeshFinalo2.println("qsub -l h_vmem=5G rt"+cntr + ".sh");
			
			cntr++;
			exesho.close();
			//exesho2.close();
			//sho.println("mkdir /data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/results/result"+cntr);
			//sho.println(backbone1+t[0]+backbone2+t[1]+backbone3+cntr++);
		}
		exeshFinalo.close();
		//exeshFinalo2.close();
		//s/ho.close();
		datain.close();
		
	}
}
