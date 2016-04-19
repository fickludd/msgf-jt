package virmid;

import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class GenerateJoinSNVMixScript {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String header = "#!/bin/bash\n#\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n#";
		String backbone1 = "jsm.py train /data/home/kyjeong/reference/Homo_sapiens_assembly19.fasta ";
		String backbone12 = "jsm.py classify /data/home/kyjeong/reference/Homo_sapiens_assembly19.fasta ";
		String backbone2 = " /data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/config/";
		String backbone3 = "--skip_size 1000 --priors_file /data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/config/priors.cfg --initial_parameters_file /data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/config/params.cfg";
		
		String data = "/home/kwj/workspace/Strelka/strelka_workflow-v0.4.10.2/scripts/ourrecal.txt";
		String trainsh = "/home/kwj/workspace/sh/st";
		String exeshFinal = "/home/kwj/workspace/sh/st.sh";
		
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
		while((s=datain.readLine())!=null){
			String[] t = s.split(" ");
			
			PrintStream exesho = new PrintStream(trainsh + cntr + ".sh");
			//PrintStream exesho2= new PrintStream(runsh + cntr + ".sh");
			
			exesho.println(header);
			exesho.println(backbone1+t[0]+" "+t[1] +backbone2 + "p"   +cntr + ".cfg "+backbone3);
			
			//exesho.println(header);
			exesho.println(backbone12+t[0]+" "+t[1]+" --model snvmix2 --out_file /data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/results" + t[0].substring(t[0].lastIndexOf("/")) + ".txt --parameters_file"+ backbone2 + "p"   +cntr + ".cfg");
			exesho.println(backbone12+t[0]+" "+t[1]+" --model snvmix2 --out_file /data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/results" + t[0].substring(t[0].lastIndexOf("/")) + ".post.txt --parameters_file"+ backbone2 + "p"   +cntr + ".cfg --somatic_threshold 0.5 --post_process");
			
			
			
			exeshFinalo.println("qsub -l h_vmem=5G -v PATH=/data/home/kyjeong/python/bin/:/data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/:$PATH -v PYTHONPATH=/data/home/kyjeong/JointSNVMix/JointSNVMix-0.8-b2/joint_snv_mix/:$PYTHONPATH st"+cntr + ".sh");
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
