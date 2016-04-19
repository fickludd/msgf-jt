package virmid;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class GenerateStrelkaScript {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String header = "#!/bin/bash\n#\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n#";
		String backbone1 = "/data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/configureStrelkaWorkflow.pl --normal=";
		String backbone2 = " --tumor=";
		String backbone3 = " --ref=/data/home/kyjeong/reference/Homo_sapiens_assembly19.fasta --config=/data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/strelka/etc/strelka_config_eland_default.ini --output-dir=/data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/results/result";
		String sh = "/home/kwj/workspace/Strelka/strelka_workflow-v0.4.10.2/scripts/conf.sh";
		String data = "/home/kwj/workspace/Strelka/strelka_workflow-v0.4.10.2/scripts/ourrecal.txt";
		String exesh = "/home/kwj/workspace/Strelka/strelka_workflow-v0.4.10.2/scripts/st";
		String exeshFinal = "/home/kwj/workspace/Strelka/strelka_workflow-v0.4.10.2/scripts/st.sh";
		
		PrintStream sho = new PrintStream(sh);
		PrintStream exeshFinalo = new PrintStream(exeshFinal);
		BufferedLineReader datain = new BufferedLineReader(data);
		String s;
		//make -C /data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/results/result1
		sho.println(header);
		int cntr = 1;
		while((s=datain.readLine())!=null){
			String[] t = s.split(" ");
			
			PrintStream exesho = new PrintStream(exesh + cntr + ".sh");
			
			exesho.println(header);
			exesho.println("make -C /data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/results/result"+cntr);
			exeshFinalo.println("qsub -l h_vmem=5G st"+cntr + ".sh");
			exesho.close();
			//sho.println("mkdir /data/home/kyjeong/Strelka/strelka_workflow-v0.4.10.2/results/result"+cntr);
			sho.println(backbone1+t[0]+backbone2+t[1]+backbone3+cntr++);
		}
		exeshFinalo.close();
		sho.close();
		datain.close();
		
	}

}
