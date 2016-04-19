package virmid;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import parser.BufferedLineReader;

public class TestSensitivity {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String[] fileNames = new String[51];
		HashMap<String, Integer> tagMap = new HashMap<String, Integer>();
		HashMap<String, HashMap<String, HashSet<Integer>>> answer = new HashMap<String, HashMap<String, HashSet<Integer>>>();
		// fil, chr, position
		int[] validated = new int[50];
		int[][] correctlyCalled = new int[5][50]; // strelka, strelkapost, jsm, jsmpost, mutect
		int[][] called = new int[5][50];
		
		String ansf = System.getProperty("user.home")+"/Dropbox/Exome/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.1.0.somatic.maf";
		String stlf= System.getProperty("user.home") + "/Dropbox/Virmid/sensitivity/StrelkaResult";
		String jsm ="/home/kwj/workspace/VirmidSensitivityResults/JointResult/results";///results";
		String mut ="/home/kwj/workspace/VirmidSensitivityResults/MetectResults_no_cosmic"; // MetectResults_no_cosmic
		String ourrecall = System.getProperty("user.home")+"/Dropbox/Exome/ourrecal.txt";
		String out = System.getProperty("user.home")+"/Dropbox/Exome/Kyowon_sensitivity";

		String s;
		
		BufferedLineReader oc = new BufferedLineReader (ourrecall);
		
		HashSet<String> files = new HashSet<String>();
		
		int i=1;
		while((s=oc.readLine())!=null){
			s = s.replaceAll("c52847b9d5bb9da3b067e9d68f3b3a01", "TCGA-E2-A14N-01A-31D-A135-09_IlluminaGA-DNASeq_exome");
			s = s.replaceAll("6bf16f633d53c28c9746de992fb66c0b", "TCGA-E2-A14N-10A-01D-A135-09_IlluminaGA-DNASeq_exome");
			
			int lastindex = s.indexOf("_Illumina");
			if(s.indexOf("_HOLD")>0){
				lastindex = lastindex > s.indexOf("_HOLD")? s.indexOf("_HOLD") : lastindex;
			}
			if(lastindex<0){
				i++;
				continue;
			}
			String k = s.substring(s.indexOf("TCGA"), lastindex);
			
			
			//System.out.println(i+ "\t" + k);
			fileNames[i] = k;
			tagMap.put(k, i-1);
			i++;
			files.add(k);
		}
		oc.close();
		
		BufferedLineReader ans = new BufferedLineReader(ansf);
		
		while((s=ans.readLine())!=null){
			if(s.startsWith("#") || s.startsWith("Hugo_Symbol")) continue;
			String[] t = s.split("\t");
			//System.out.println(t[9] + "\t" + t[24] + "\t" + t[15]);
			if(!t[9].equals("SNP") || !t[24].equals("Valid") || !files.contains(t[15])) continue;
			
			if(!answer.containsKey(t[15])) answer.put(t[15], new HashMap<String, HashSet<Integer>>());
			if(!answer.get(t[15]).containsKey(t[4])) answer.get(t[15]).put(t[4], new HashSet<Integer>());
			int po = Integer.parseInt(t[5]);
			//if(t[2].equals("MT"))System.out.println(t[2]);
			answer.get(t[15]).get(t[4]).add(po);
			validated[tagMap.get(t[15])]++;
		}
		ans.close();

		//for(int h=0;h<50;h++)
		//	System.out.println(validated[h]);
		
	
		
		File stdir = new File(stlf);
		File jsmdir = new File(jsm);
		File mutdir = new File(mut);
		
		//PrintStream outs = new PrintStream(out);
		int thresholdStrelka = 0;
		PrintStream tb = new PrintStream(out+ thresholdStrelka + "_DP.table.txt");
		PrintStream answeredStrelka = new PrintStream(out + "_StrelkaAnswer.txt");
		PrintStream answeredStrelka2 = new PrintStream(out + "_StrelkaAnswer_filtered.txt");
		
		PrintStream answeredMutect = new PrintStream(out + "_MutectAnswer.txt");
		PrintStream answeredJSM = new PrintStream(out + "_JSMAnswer.txt");
		PrintStream answeredJSM2 = new PrintStream(out + "_JSMAnswer_filtered.txt");
		
		//outs.println("#No\tFile\tChr\tPos");
		//outs2.println("#No\tFile\tChr\tPos");
		
		tb.println("TAG\t# validated SNP\t# Strelka called SNP\t# Strelka correctly called SNP\tStrelka Sensitivity (%)" +
				"\t# Strelka (filtered) called SNP\t# Strelka (filtered) correctly called SNP\tStrelka (filtered) Sensitivity (%)" +
				"\t# JSM called SNP\t# JSM correctly called SNP\tJSM Sensitivity (%)" +
				"\t# JSM (filtered) called SNP\t# JSM (filtered) correctly called SNP\tJSM (filtered) Sensitivity (%)" +
				"\t# MuTect called SNP\t# MuTect correctly called SNP\tMuTect Sensitivity (%)");
		//tb2.println("TAG\t# validated SNP\t# called SNP\t# correctly called SNP\tSensitivity (%)");
		
		for(File stf : stdir.listFiles()){
			if(stf.isFile()) continue;
			String path = stf.getAbsolutePath();
			String k = fileNames[Integer.parseInt(path.substring(path.lastIndexOf('t')+1))];
			//System.out.println(Integer.parseInt(path.substring(path.lastIndexOf('t')+1)) + "\t" + k);
			HashMap<String, HashSet<Integer>> answer1 = answer.get(k);
			if(answer1 == null){
				//System.out.println("No.. "+k);
				continue;
			}
			HashMap<String, HashMap<Integer, Boolean>> answered = new HashMap<String, HashMap<Integer, Boolean>>();
			
			for(String ss : answer1.keySet()){
				answered.put(ss, new HashMap<Integer, Boolean>());
				for(int iii : answer1.get(ss)){
					answered.get(ss).put(iii, false);
				}
			}
			
			//##FILTER=<ID=QSS_ref,Description="Normal sample is not homozygous ref or ssnv Q-score < 15, ie calls with NT!=ref or QSS_NT < 15">
			if(new File(path + "/results/all.somatic.snvs.vcf").exists()){
				BufferedLineReader sr = new BufferedLineReader(path + "/results/all.somatic.snvs.vcf");
				while((s=sr.readLine())!=null){
					if(s.startsWith("#"))continue;
					String[] t = s.split("\t");
					
					boolean pass = true;
					if(t[6].equals("PASS") || t[6].equals("DP")) pass = true;
					else pass =false;
					/*if(t[6].equals("QSS_ref") || t[6].equals("QSS_ref;DP")){
						String[] tt = t[7].split(";");
						for(String filter : tt){
							if(filter.startsWith("QSS_NT")){
								
								int score = Integer.parseInt(filter.substring(filter.lastIndexOf('=')+1));
								//System.out.println(filter + "\t" + score);
								if(score < thresholdStrelka){
									pass = false;
									break;
								}
							}
						}
					}else if(t[6].equals("DP"))pass = true;
					else if(t[6].startsWith("PASS")) pass =true;
					else pass = false;
					*/
					if(!pass) continue;
					
					HashSet<Integer> answer2 = answer1.get(t[0]);
					if(answer2 != null){
						//answered.get(t[0]).put(Integer.parseInt(t[1]), true);
						if(answer2.contains(Integer.parseInt(t[1]))){
							answeredStrelka.println(k+"\t"+t[0]+"\t"+t[1]);
							//outs.println(correctlyCalled[0][tagMap.get(k)]+1 + "\t" + k + "\t" + t[0] + "\t" + t[1]);						
							correctlyCalled[0][tagMap.get(k)]++;
						}
					}
					called[0][tagMap.get(k)]++;
				}
				tb.print(k+"\t"+validated[tagMap.get(k)]);
				//System.out.println(k);
				//outs.println("\n  ## Validated SNP: " + validated[tagMap.get(k)] + "\tCalled SNP: " + called[0][tagMap.get(k)] + "\tCorrectly called SNP: " + correctlyCalled[0][tagMap.get(k)]+"\n");
				
				sr.close();
			}
			if(new File(path + "/results/passed.somatic.snvs.vcf").exists()){
				BufferedLineReader srf = new BufferedLineReader(path + "/results/passed.somatic.snvs.vcf");
				while((s=srf.readLine())!=null){
					if(s.startsWith("#"))continue;
					String[] t = s.split("\t");
					HashSet<Integer> answer2 = answer1.get(t[0]);
					if(answer2 != null){
						if(answer2.contains(Integer.parseInt(t[1]))){
							answeredStrelka2.println(k+"\t"+t[0]+"\t"+t[1]);
							
							correctlyCalled[1][tagMap.get(k)]++;
						//	outs.println(correctlyCalled[1][tagMap.get(k)]+1 + "\t" + k + "\t" + t[0] + "\t" + t[1]);						
							
						}
					}
					called[1][tagMap.get(k)]++;
				}
				
			//	outs2.println("\n  ## Validated SNP: " + validated[tagMap.get(k)] + "\tCalled SNP: " + called[1][tagMap.get(k)] + "\tCorrectly called SNP: " + correctlyCalled[1][tagMap.get(k)]+"\n");
				
				
				srf.close();
			}
			for(File jsmf : jsmdir.listFiles()){

				if(jsmf.isDirectory()) continue;
				String fn = jsmf.getAbsolutePath();
				if(fn.endsWith("post.txt"))continue;
				String fn1 = fn.replaceAll("c52847b9d5bb9da3b067e9d68f3b3a01", "TCGA-E2-A14N-01A-31D-A135-09_IlluminaGA-DNASeq_exome");
				String fn2 = fn1.replaceAll("6bf16f633d53c28c9746de992fb66c0b", "TCGA-E2-A14N-10A-01D-A135-09_IlluminaGA-DNASeq_exome");
				if(fn2.contains(k)){
					//System.out.println(fn2);
					BufferedLineReader in = new BufferedLineReader(fn);
					String a;
					
					while((a=in.readLine())!=null){
						if(a.startsWith("chrom"))continue;
						String[] t = a.split("\t");
						//String ch = tok[0];
						float p=Float.parseFloat(t[9])+Float.parseFloat(t[10]);
						HashSet<Integer> answer2 = answer1.get(t[0]);
						//if(answer2 != null && answer2.contains(Integer.parseInt(t[1]))){
						//	System.out.println(p+"\t"+a);
						//}
						
						if(p>=.5f){
							//HashSet<Integer> answer2 = answer1.get(t[0]);
							if(answer2 != null){
								answered.get(t[0]).put(Integer.parseInt(t[1]), true);
								if(answer2.contains(Integer.parseInt(t[1]))){
									answeredJSM.println(k+"\t"+t[0]+"\t"+t[1]);
									
									correctlyCalled[2][tagMap.get(k)]++;
								}
							}
							called[2][tagMap.get(k)]++;
						}
					}
					
					in.close();
				}
				
			}
			
			
			for(File jsmf : jsmdir.listFiles()){

				if(jsmf.isDirectory()) continue;
				String fn = jsmf.getAbsolutePath();
				if(!fn.endsWith("post.txt"))continue;
				String fn1 = fn.replaceAll("c52847b9d5bb9da3b067e9d68f3b3a01", "TCGA-E2-A14N-01A-31D-A135-09_IlluminaGA-DNASeq_exome");
				String fn2 = fn1.replaceAll("6bf16f633d53c28c9746de992fb66c0b", "TCGA-E2-A14N-10A-01D-A135-09_IlluminaGA-DNASeq_exome");
				if(fn2.contains(k)){
					//System.out.println(fn2);
					BufferedLineReader in = new BufferedLineReader(fn);
					String a;
					
					while((a=in.readLine())!=null){
						if(a.startsWith("chrom"))continue;
						String[] t = a.split("\t");
						//String ch = tok[0];
						float p=Float.parseFloat(t[9])+Float.parseFloat(t[10]);
						if(p>=.5f){
							HashSet<Integer> answer2 = answer1.get(t[0]);
							if(answer2 != null){
								//answered.get(t[0]).put(Integer.parseInt(t[1]), true);
								if(answer2.contains(Integer.parseInt(t[1]))){
									answeredJSM2.println(k+"\t"+t[0]+"\t"+t[1]);
									
									correctlyCalled[3][tagMap.get(k)]++;
								}
							}
							called[3][tagMap.get(k)]++;
						}
					}
					
					in.close();
				}
			}
			
			
			
			for(File mutf : mutdir.listFiles()){
				if(mutf.isDirectory()) continue;
				String fn = mutf.getAbsolutePath();
				if(fn.endsWith("cov.txt"))continue;
				String fn1 = fn.replaceAll("c52847b9d5bb9da3b067e9d68f3b3a01", "TCGA-E2-A14N-01A-31D-A135-09_IlluminaGA-DNASeq_exome");
				String fn2 = fn1.replaceAll("6bf16f633d53c28c9746de992fb66c0b", "TCGA-E2-A14N-10A-01D-A135-09_IlluminaGA-DNASeq_exome");
				if(fn2.contains(k)){
					//System.out.println(fn2);
					BufferedLineReader in = new BufferedLineReader(fn);
					String a;
					
					while((a=in.readLine())!=null){
						if(a.startsWith("##")||a.startsWith("contig"))continue;
						String[] t = a.split("\t");
						//String ch = tok[0];
						
						if(!a.endsWith("REJECT")){
							HashSet<Integer> answer2 = answer1.get(t[0]);
							//if(!answered.containsKey(t[0])) answered.put(t[0], new HashMap<Integer, Boolean>());

							if(answer2 != null){
								//System.out.println(a + "\t" + fn);
								//answered.get(t[0]).put(Integer.parseInt(t[1]), true);
									
								
								if(answer2.contains(Integer.parseInt(t[1]))){
									answeredMutect.println(k+"\t"+t[0]+"\t"+t[1]);
									
									correctlyCalled[4][tagMap.get(k)]++;
								}
							}
							called[4][tagMap.get(k)]++;
						}
					}
					
					in.close();
				}
			}
			
			for(String key : answered.keySet()){
				for(int keyi : answered.get(key).keySet()){
					if(!answered.get(key).get(keyi)){
						System.out.println( k+"\t"+key+"\t"+keyi);
					}
				}
			}
			
			for(int j=0;j<5;j++)
				tb.print("\t"+called[j][tagMap.get(k)]+"\t"+ correctlyCalled[j][tagMap.get(k)]+"\t"+String.format("%.1f", (float)100f*correctlyCalled[j][tagMap.get(k)]/validated[tagMap.get(k)]));
			tb.println();
		}
		
		
		
		//System.out.println();
		
		int[] t = new int[3];
		for(int h=0;h<50;h++)
			t[0] += validated[h];
		
		for(int h=0;h<50;h++)
			t[1] += called[0][h];
		
		for(int h=0;h<50;h++)
			t[2] += correctlyCalled[0][h];
		
		//outs.println("\n  ### Total validated SNP: " + t[0] + "\tTotal called SNP: " + t[1] + "\tTotal Correctly called SNP: " + t[2]+"\n");
		tb.print("Total\t"+t[0]+"\t"+t[1]+"\t"+t[2]+"\t"+String.format("%.1f", (float)100f*t[2]/t[0]));
		t = new int[3];
		for(int h=0;h<50;h++)
			t[0] += validated[h];
		
		for(int h=0;h<50;h++)
			t[1] += called[1][h];
		
		for(int h=0;h<50;h++)
			t[2] += correctlyCalled[1][h];
		tb.println("\t"+t[1]+"\t"+t[2]+"\t"+String.format("%.1f", (float)100f*t[2]/t[0]));
		//outs2.println("\n  ### Total validated SNP: " + t[0] + "\tTotal called SNP: " + t[1] + "\tTotal Correctly called SNP: " + t[2]+"\n");
		
		tb.close();
		answeredStrelka.close();
		answeredMutect.close();
		answeredStrelka2.close();
		answeredJSM.close();
		answeredJSM2.close();
		
		//tb2.close();
		//outs2.close();
		
		//outs.close();
		
		
	}

}
