package uninovoOld.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.Composition;
import msutil.Peptide;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import uninovoOld.DenovoReconstruction;
import uninovoOld.Node;

public class DrawVennDiagram {

	
	private static float[] getUniNovoResults(HashSet<String> ret, String file, String specFile, int charge, int numRec) throws IOException{
		BufferedLineReader in = new BufferedLineReader(file);
		
		String s;
		float len = 0;
		float pepLen = 0;
		float tp =0, fp = 0;
		
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(spec.getCharge() != charge) continue;
			pepLen += spec.getAnnotationStr().length()-1;
		}
		
		int qspec = 0;
		int cspec = 0;
		while((s = in.readLine())!=null){
			if(s.startsWith("##")){
				String[] t = s.split("\t");
				if(t.length > 4){
					tp = Float.parseFloat(t[1]);
					fp = Float.parseFloat(t[2]);
					
					qspec = Integer.parseInt(t[3]);
					cspec = Integer.parseInt(t[4]);
				}
				continue;
			}
			if(Integer.parseInt(s.split("\t")[0]) != charge) continue;
			
			String[] t = s.split("\t");
			String k = t[0] + "\t" + t[1];

			len += Float.parseFloat(t[2]);
			while(ret.contains(k)){
				k+="'";
			}
			
			ret.add(k);
			
			
		}	
		
		in.close();
		
		float[] f = new float[2];
		f[0] = len/ret.size();
		f[1] = qspec;
		
		if(numRec==1){
			float pre =  (tp/(tp+fp));
			float rec =  (tp/pepLen);
			float fscore = (2*pre*rec)/(pre+rec);
			System.out.println("%UniNovo Precision, Recall, Fscore : \n" + pre  + " " + rec + " " + fscore);
		}
		
		return f;
	}
	
	private static float getPepNovoResults(HashSet<String> ret, String file, String specFile, int charge, int numRec) throws IOException{
		BufferedLineReader in = new BufferedLineReader(file);
		String s;
		// = new  HashSet<String>();
	
		String annotation = null;
		boolean newSpec = false, start = false, isCorrect = false;
		float len = 0;
		float tp = 0;
		float fp = 0;
		float pepLen = 0;
		
		ArrayList<Float> prms = null;
		int cutlen = 0;
		ArrayList<String> recsLIQK = null;
		
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(spec.getCharge() != charge) continue;
			
			pepLen += spec.getAnnotationStr().length()-1;
		}
		
		while((s = in.readLine())!=null){
			if(s.startsWith("#") || s.isEmpty()) continue;
			if(s.startsWith(">>")){
				annotation = s.split(" ")[s.split(" ").length-3];
				newSpec = true; start = true;isCorrect = false;
				recsLIQK = new ArrayList<String>();
				if(prms != null){
					cutlen -= prms.size();
					tp += cutlen;
				}
				prms = new ArrayList<Float>();
			//	System.out.println(annotation);
				for(float prm : new Peptide(annotation).getPRMMasses(true, 0))
					prms.add(prm);
				cutlen = prms.size();
				
				continue;
			}
			
			if(!start) continue;
			
			String[] t = s.split("\t");
		//	System.out.println(s);
			if(charge > 0 && Integer.parseInt(t[7]) != charge) continue;
			
			if(newSpec){
				newSpec = false;
			}
			
			String p = t[8].replace('L', 'I');
			
			if(tol.getToleranceAsDa(100) > 0.3f)
				p = p.replace('Q', 'K');
			
			if(!recsLIQK.contains(p))
				recsLIQK.add(p);
			
			if(recsLIQK.size() > numRec || isCorrect) continue;
			float cGap =  Float.parseFloat(t[5]);
			float nGap = Float.parseFloat(t[4]);
			if(cGap != 0f) cGap  -=  Composition.H2O + Composition.PROTON;//  19.017f;
			DenovoReconstruction r = new DenovoReconstruction(nGap , cGap , new Peptide(t[8]),  tol,  pmtol);
			for(int m=0; m <r.getNodes().size(); m++){
				Node node = r.getNodes().get(m);
				if(node.getMass() == 0 || node.isSink()) continue;
				
				boolean c = false;
				
				int i = Collections.binarySearch(prms, node.getMass());
				if(i>=0){
					prms.remove(i);
					c = true;
				}
				else{
					i=-i-1;
					
					for(int j=-1;j<=0;j++){
						if(i+j>=0 && i+j < prms.size())
							if(Math.abs(prms.get(i+j) - node.getMass()) < node.getTol().getToleranceAsDa(Math.max(node.getMass(), r.getPeptideMass() -node.getMass()))){
								c = true;
								prms.remove(i+j);
								break;
							}
					}
				}
				if(!c) fp ++;
			}
			
	//		if(!r.isCorrect(new Peptide(annotation)) && s.endsWith("*")){
	//			System.out.println("** " + s + "\t" + annotation);
		//	}
			//r.isCorrect(new Peptide(annotation)) && 
		//	if(r.isCorrect(new Peptide(annotation)) && !s.endsWith("*")){
		//			System.out.println("** " + s + "\t" + annotation);
		//	}
						
				
			if(!isCorrect &&  r.isCorrect(new Peptide(annotation))){//s.endsWith("*")){
				isCorrect = true;
			//	System.out.println(annotation);
				String k = charge+"\t" +annotation;
				len += t[8].length() + (Float.parseFloat(t[4]) > 0? 1 : 0) + (Float.parseFloat(t[5]) > 0? 1 : 0);

				while(ret.contains(k)){
					k+="'";
				}
				
				ret.add(k);
				
			}
		}
		in.close();
		if(numRec==1){
			float pre =  (tp/(tp+fp));
			float rec =  (tp/pepLen);
			float fscore = (2*pre*rec)/(pre+rec);
			System.out.println("%PepNovo Precision, Recall, Fscore : \n" + pre  + " " + rec + " " + fscore);
		}
	//	System.out.println("PepNovo+ : " + ret.size() + "\t" + len/ret.size() + "\t" + correctCutLen + "\t" + incorrectCutLen + "\t" + (correctCutLen)/(correctCutLen + incorrectCutLen));
		return len/ret.size();
	}
	
	private static float getPeaksResults(HashSet<String> ret, String file, String specFile, int charge, int numRec, float threshold) throws IOException{
		//HashSet<String> ret = new  HashSet<String>();
		
		BufferedLineReader in = new BufferedLineReader(file);
		String s;
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		
		HashMap<String, String> annotationMap = new HashMap<String, String>();
		
		int ll = 0;
		float pepLen = 0;
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(spec.getCharge() != charge) continue;
			pepLen += spec.getAnnotationStr().length()-1;
			
			String key = String.format("%.4f", spec.getPrecursorPeak().getMz());
			//if(spec.getAnnotationStr().equals("HEQISDLER")) System.out.println(key + "\t" + (spec.getPrecursorPeak().getMz()-(5e-5f)));
			annotationMap.put(key, spec.getAnnotationStr());
	
		}

		ArrayList<DenovoReconstruction> recs = null;
		
		int prev = -1, sn = 0;
		float len = 0;
		float tp = 0, fp = 0;
		boolean isCorrect = false;
		
		String key = null;
		while((s=in.readLine())!=null){
		
			if(s.startsWith("Scan")) continue;
			String[] t = s.split(",");
			if(t.length<6) continue;
			//System.out.println(t[5]+"\t" + charge);
			if(Integer.parseInt(t[5])!=charge) continue;
		
			int scanNum = Integer.parseInt(t[0].split(":")[t[0].split(":").length-1]);

			
			if(prev != scanNum){ // init
				isCorrect = false;
				if(recs != null){
					if(recs.size() < numRec) ll++;
					String annotation = annotationMap.get(key);
					if(annotation == null);// System.out.println(key);
					else{
						
						ArrayList<Float> prms = new ArrayList<Float>();
						for(float prm : new Peptide(annotation).getPRMMasses(true, 0))
							prms.add(prm);
						int cutlen = prms.size();
							
						for(DenovoReconstruction r : recs){
		
							for(Node node : r.getNodes()){
								if(node.getMass() == 0 || node.isSink()) continue;
								boolean c = false;
								
								int i = Collections.binarySearch(prms, node.getMass());
								if(i>=0){
									prms.remove(i);
									c = true;
								}
								else{
									i=-i-1;
									
									for(int j=-1;j<=0;j++){
										if(i+j>=0 && i+j < prms.size())
											if(Math.abs(prms.get(i+j) - node.getMass()) < node.getTol().getToleranceAsDa(Math.max(node.getMass(), r.getPeptideMass() -node.getMass()))){
												prms.remove(i+j);
												c = true;
												break;
											}
									}
								}
								if(!c) fp++;
							}

							//if(r.isCorrect(new Peptide(annotation))){
							//if(rec.equals(annotation.replace('K', 'Q').replace('L', 'I'))){
							//System.out.println(r + "\t" + annotation + "\t" + new Peptide(annotation).getMass());
							if(!isCorrect && r.isCorrect(new Peptide(annotation))){
								isCorrect = true;
								String k = charge+"\t" +annotation;
								len += r.length();
								while(ret.contains(k)){
									k+="'";
								}
								
								ret.add(k);
								break;
							}
						}
						
						cutlen -= prms.size();
						tp += cutlen;
					}
				}

				key = t[4];
				recs = new ArrayList<DenovoReconstruction>();
	
				prev = scanNum;
				sn++;
			}
			
			if(recs.size() >= numRec) continue;
			
			Peptide pep = new Peptide(t[1].replace("C(+57.02)", "C"));
			ArrayList<Float> l = new ArrayList<Float>();
			
			for(String k : t[t.length-1].split(" ")){
				l.add(Float.parseFloat(k));
			}
		//	System.out.println(pep + "\n" + l);
		//	System.exit(0);
			
			boolean isPrevGap = false;
			ArrayList<Node> nodes = new ArrayList<Node>();
			nodes.add(new Node(0f, tol, pep.getMass(), 0));
			float sum = 0;
			for(int k=0; k<pep.size()-1; k++){
				float score = l.get(k);
				if(score < threshold){
					isPrevGap = true;
				}
				else{
					if(isPrevGap){
						nodes.add(new Node(sum, tol, pep.getMass(), 0));
					}
					nodes.add(new Node(sum+pep.get(k).getMass(), tol, pep.getMass(), 0));
					isPrevGap = false;
				}
				
				sum += pep.get(k).getMass();
				
			}
			if(l.get(l.size()-1) >= threshold && isPrevGap){
				nodes.add(new Node(sum, tol, pep.getMass(), 0));
			}
			
			Node sink = new Node(pep.getMass(),  pmtol, pep.getMass(), 0);
			sink.setSink(pmtol);
			nodes.add(sink);
			DenovoReconstruction rec = new DenovoReconstruction(nodes);
			if(!recs.contains(rec)) recs.add(rec);
			//System.out.println(pep + "\t" + pep.getMass() + "\n" + l + "\n" + new DenovoReconstruction(nodes) + "\t" + nodes.get(nodes.size()-1).getMass() + "\t" + new DenovoReconstruction(nodes) .length());
		//	System.exit(0);
		}
	
		
		in.close();
		
		if(numRec==1){
			float pre =  (tp/(tp+fp));
			float rec =  (tp/pepLen);
			float fscore = (2*pre*rec)/(pre+rec);
			System.out.println("%PEAKS Precision, Recall, Fscore : \n" + pre  + " " + rec + " " + fscore);
		}
		
		//System.out.println("PEAKS : " + ret.size() + "\t" + len/ret.size() + "\t" + correctCutLen + "\t" + incorrectCutLen + "\t" + (correctCutLen)/(correctCutLen + incorrectCutLen));
		if(ll>0) System.out.println("%" + ll);
		return len/ret.size() ;
	}
	
	private static float getpNovoResults(HashSet<String> ret, String file, String specFile, int charge, int numRec) throws IOException{
		//HashSet<String> ret = new  HashSet<String>();
		
		BufferedLineReader in = new BufferedLineReader(file);
		String s;
		Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
		
		HashMap<String, String> annotationMap = new HashMap<String, String>();
		
		float pepLen = 0;
		while(iterator.hasNext()){
			Spectrum spec = iterator.next();
			if(spec.getCharge() != charge) continue;
			pepLen += spec.getAnnotationStr().length()-1;
			
			String key = spec.getTitle();
		//	if(spec.getAnnotationStr().equals("LFLEPTRK")) System.out.println(key);
			annotationMap.put(key, spec.getAnnotationStr());
	
		}
		
		ArrayList<String> recs = null;
		ArrayList<String> recsLIQK = null;
		
		int sn = 0;
		float len = 0;
		float tp = 0, fp = 0;
		int cn = -1;
		boolean isCorrect = false;
		
		String key = null;
		while((s=in.readLine())!=null){
			if(s.isEmpty()) continue;
			String[] t = s.split("\t");
			if(cn >0 && cn != charge) continue;
			if(s.startsWith("Elution")){ // init
				isCorrect = false;
				cn = Integer.parseInt(t[2]);
				if(cn != charge){
					cn = -1;
				}
				if(recs != null){
					String annotation = annotationMap.get(key);
					//System.out.println(recs.size());
					if(annotation == null);// System.out.println(key);
					else{
						ArrayList<Float> prms = new ArrayList<Float>();
						for(float prm : new Peptide(annotation).getPRMMasses(true, 0))
							prms.add(prm);
						int cutlen = prms.size();
						
						for(String rec : recs){
							DenovoReconstruction r = new DenovoReconstruction(0, 0, new Peptide(rec),  tol,  pmtol);
							for(Node node : r.getNodes()){
								if(node.getMass() == 0 || node.isSink()) continue;
								boolean c = false;
								
								int i = Collections.binarySearch(prms, node.getMass());
								if(i>=0){
									prms.remove(i);
									c = true;
								}
								else{
									i=-i-1;
									
									for(int j=-1;j<=0;j++){
										if(i+j>=0 && i+j < prms.size())
											if(Math.abs(prms.get(i+j) - node.getMass()) < node.getTol().getToleranceAsDa(Math.max(node.getMass(), r.getPeptideMass() -node.getMass()))){
												prms.remove(i+j);
												c = true;
												break;
											}
									}
								}
								if(!c) fp++;
							}

							//if(r.isCorrect(new Peptide(annotation))){
							//if(rec.equals(annotation.replace('K', 'Q').replace('L', 'I'))){
							if(!isCorrect && r.isCorrect(new Peptide(annotation))){
								if(!annotation.toString().replace('L', 'I').equals(rec)){
									break;
								}
								isCorrect = true;
								String k = charge+"\t" +annotation;
								len += rec.length();
								while(ret.contains(k)){
									k+="'";
								}
								
								ret.add(k);
								break;
							}
						}
						
						cutlen -= prms.size();
						tp += cutlen;
					}
				}

				key = t[0];//String.format("%.4f", Float.parseFloat(t[1]));
				recs = new ArrayList<String>();
				recsLIQK = new ArrayList<String>();
				sn++;
				continue;
			}
			
			if(recsLIQK.size() >= numRec) continue;
			String rec = t[1].replace('L', 'I');
			rec = rec.substring(2,rec.length()-2);
			//System.out.println(rec);
			if(!recs.contains(rec)){
				recs.add(rec);
			}
			String recLIQK = rec;
			if(tol.getToleranceAsDa(100) > 0.3f)
				recLIQK = rec.replace('Q', 'K');
			if(!recsLIQK.contains(recLIQK)){
				recsLIQK.add(recLIQK);
			}
		}
		
		in.close();
		if(numRec==1){
			float pre =  (tp/(tp+fp));
			float rec =  (tp/pepLen);
			float fscore = (2*pre*rec)/(pre+rec);
			System.out.println("%PNovo Precision, Recall, Fscore : \n" + pre  + " " + rec + " " + fscore);
		}
		
		return len/ret.size();
	}

	private static int[] getOverlaps(HashSet<String> a, HashSet<String> b, HashSet<String> c){
		int[] i = new int[7];
		
		for(String s : a){
			if(b.contains(s)){ // ab
				i[3]++;
				if(c.contains(s)) // abc 
					i[6]++;
			}
			if(c.contains(s)){ // ca
				i[5]++;
			}
		}
		
		for(String s : b){ 
			if(c.contains(s)){ // bc
				i[4]++;
			}
		}
		
		i[3] = i[3] - i[6]; //ab
		i[4] = i[4] - i[6]; //bc
		i[5] = i[5] - i[6]; //ca
		
		i[0] = a.size() - i[3]-i[5]-i[6];
		i[1] = b.size() - i[3]-i[4]-i[6];
		i[2] = c.size() - i[4]-i[5]-i[6];
		
		return i;
	}
	
	
	static Tolerance tol = new Tolerance(0.5f, false);
	static Tolerance pmtol = new Tolerance(0.5f, false);

	public static void main(String[] args) throws IOException {
		int charge = 2;
		
		float peaksThreshold = 30;
		String frag = "STD";
		String enzyme = "Trypsin";//"";
		float setAccuracy = 0.0f;
		boolean paired = false;
		
		int[] numRecs = {1,5,20};
		
		int[][] crecs = new int[5][25];
		float[][] clen = new float[5][25];
		
		System.out.println(frag+enzyme+charge+"=[");
		
		for(int numRec : numRecs){
			if(frag.equals("HCD")) enzyme = "";
			
			String specFile = "/home/kwj/Dropbox/Test/" + frag  + enzyme +"_test.mgf";//Trypsin_test
			//String specFile = "/home/kwj/Dropbox/Test/CIDETDTrypsin/CID.mgf";
			String pepNovo =  "/home/kwj/Dropbox/"+ frag + enzyme + "PepNovo.txt";
			String UniNovo = "/home/kwj/Dropbox/Test/" + frag  + enzyme + "_test.mgf" + (paired? "paired": "" )+ ".out" + numRec +"_"+ charge;
			String peaks = "/home/kwj/Dropbox/"+frag+enzyme+"_DENOVO/all de novo candidates.csv";
			String MSGFDB = "/home/kwj/Dropbox/Test/" + frag + enzyme + "_testPRMMSGF.mgf.out" + numRec +"_"+ charge;
			String pNovo = "/home/kwj/Dropbox/Test/HCD_pNovo.txtPNOVO_TEMP.0";
			
			if(frag.equals("STD")){
				specFile =  "/home/kwj/Dropbox/Test/Standard1388_corrected.mgf";
				pepNovo = "/home/kwj/Dropbox/Test/Standard.txt";
				peaks = "/home/kwj/Dropbox/Standard2_DENOVO_2/all de novo candidates.csv";
				UniNovo = "/home/kwj/Dropbox/Test/Standard1388_corrected.mgf.out" + numRec +"_"+ charge;
			}
			
			
			if(setAccuracy > 0) UniNovo = UniNovo + "_" + setAccuracy;
			
			//String UniNovo = "/home/kwj/Dropbox/Test/CIDETDTrypsin.out" + numRec +"_"+ charge;
			
			ArrayList<String> all = new ArrayList<String>();
			
			ArrayList<HashSet<String>> deNovoResults = new ArrayList<HashSet<String>>();
				
			HashMap<BitSet, Integer> venn = new HashMap<BitSet, Integer>();
	
			if(frag.equals("HCD"))  tol = new Tolerance(20f, true);
			
			HashSet<String> uret = new HashSet<String>();
			float[] ss = getUniNovoResults(uret, UniNovo, specFile, charge, numRec);
			clen[0][numRec] = ss[0];
			crecs[0][numRec] = uret.size();
			deNovoResults.add(uret);
			
			//deNovoResults.add();
			
			HashSet<String> peaksret = new HashSet<String>();
			HashSet<String> pepnovosret = new HashSet<String>();
			HashSet<String> pnovosret = new HashSet<String>();
			
			if(setAccuracy == 0){
				
				if(!frag.equals("CIDETD") && !frag.equals("HCD")){
					
					clen[2][numRec] = getPeaksResults(peaksret, peaks, specFile, charge, numRec,peaksThreshold);
					crecs[2][numRec] = peaksret.size();
					deNovoResults.add(peaksret);
				}
				if(!frag.equals("ETD")){
					clen[1][numRec] = getPepNovoResults(pepnovosret, pepNovo, specFile, charge, numRec);
					crecs[1][numRec] = pepnovosret.size();
					deNovoResults.add(pepnovosret);
				}
				if(frag.equals("HCD")){
					clen[3][numRec] = getpNovoResults(pnovosret, pNovo, specFile, charge, numRec);
					crecs[3][numRec] = pnovosret.size();
					deNovoResults.add(pnovosret);
				}
				
//				if(!frag.equals("STD")){
//					uret = new HashSet<String>();
//					clen[4][numRec] = getUniNovoResults(uret, MSGFDB, specFile, charge, numRec)[0];
//					crecs[4][numRec] = uret.size();
//					deNovoResults.add(uret);
//				}
			}
			
			if(numRec == 1) System.out.println("];");
			
			Iterator<Spectrum> iterator = new SpectraIterator(specFile, new MgfSpectrumParser());
			
			HashSet<String> tmp = new HashSet<String>();
			
			while(iterator.hasNext()){
				Spectrum s = iterator.next();
				if(s.getCharge() != charge) continue;
				String t = charge+ "\t" + s.getAnnotationStr();
				
				while(tmp.contains(t)){
					t+="'";
				}
				
				all.add(t);
				tmp.add(t);
				
			}
			//System.out.println("Total spec Num : " + all.size());
			
			
			
			for(String rec : all){
				BitSet key = new BitSet();
				for(int i=0; i< deNovoResults.size(); i++){
					HashSet<String> recs = deNovoResults.get(i);
					if(recs.contains(rec)) key.set(i);	
				}
				
				if(!venn.containsKey(key)) venn.put(key, 0);
				venn.put(key, venn.get(key)+1);
			}
			
			if(setAccuracy > 0){
			//	System.out.println("q"  + frag + (paired? "paired": "" ) + enzyme + charge +"_" + numRec + "=" + ss[1] + ";");
			}
			
			//ss
			
			int ov[] = getOverlaps(uret, pepnovosret, pnovosret);
			
			if(!frag.equals("HCD")){
				ov = getOverlaps(uret, pepnovosret, peaksret);
			}
			
			for(int o : ov) System.out.print(o+"\t");
			System.out.println();
			System.out.print("n" + frag + (paired? "paired": "" ) + enzyme + charge +"_"+ numRec + "=[");
			for(int i=0;i<crecs.length;i++){
				System.out.print(crecs[i][numRec]+" ");
			}
			System.out.println("];");
			

			
			
			
		}
		
		
		
		//HashSet<String> peaksret = new HashSet<String>();
		//HashSet<String> pepnovosret = new HashSet<String>();
		//HashSet<String> pnovosret
		
		
		
		for(int numRec : numRecs){
			System.out.print("l"  + frag + (paired? "paired": "" ) + enzyme + charge +"_" + numRec + "=[");
			for(int i=0;i<crecs.length;i++){
				System.out.print(clen[i][numRec]+" ");
			}
			System.out.println("];");
		}
		
	}

}
