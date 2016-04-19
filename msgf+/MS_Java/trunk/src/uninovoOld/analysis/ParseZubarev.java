package uninovoOld.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

import msgf.NominalMass;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import parser.BufferedLineReader;

public class ParseZubarev {
        
        static AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
        
        static private class reconstructions{
                String title;
                //int sn;
                ArrayList<Integer> modes;
                HashMap<Integer, ArrayList<ArrayList<Integer>>> recs;
                HashMap<Integer, ArrayList<BitSet>> isGap;
                
                public boolean hit = false;
                public boolean PTMhit = false;
                
                reconstructions(String title, int mode, ArrayList<String> recs){
                        this.title = title;
                        
                        this.modes = new ArrayList<Integer>();
                        this.modes.add(mode);
                        this.recs = new HashMap<Integer, ArrayList<ArrayList<Integer>>>();
                        this.isGap = new HashMap<Integer, ArrayList<BitSet>>();
                        
                        ArrayList<ArrayList<Integer>> t = new ArrayList<ArrayList<Integer>>();
                        ArrayList<BitSet> gs = new  ArrayList<BitSet>();
                        for(String r : recs){
                        	BitSet g = new BitSet();
                        	t.add(getrecfrom(r, g));
                        	gs.add(g);
                        }
                        this.recs.put(mode, t);
                        this.isGap.put(mode, gs);
                }
                
                private ArrayList<Integer> getrecfrom(String rec, BitSet isGap){
                        ArrayList<Integer> out = new ArrayList<Integer>();
                        ArrayList<Integer> tmp = new ArrayList<Integer>();
                        int m = 0;
                        String t = "";
                        
                        int j = 0;
                        for(int i=0; i<rec.length(); i++){
                            char c = rec.charAt(i);
                            if(c == '['){
                                    m = 1;
                                    t = "";
                                    continue;
                            }else if(c == ']'){
                                    m = 0;
                              //      if(!tmp.isEmpty()) //499	3264	40368
                                    	tmp.add(NominalMass.toNominalMass(Float.parseFloat(t)));
                                    	isGap.set(j);
                                    	j++;
                                    continue;
                            }else if(m == 1){
                                    t += c;
                            }else{
                                    tmp.add(aaSet.getAminoAcid(c).getNominalMass());
                                    j++;
                            }
                                
                        }
                        
                        int pre = 0;
                        for(int x : tmp){
                                out.add(pre + x);
                                pre += x;
                        }
                     //   System.out.println(rec + "\t" + out);
                        return out;
                }
                
                void add(reconstructions r){
                	for(int mode : r.modes){                		
                		if(!this.recs.containsKey(mode)){
                			this.modes.add(mode);
                			this.recs.put(mode, new ArrayList<ArrayList<Integer>>());
                			this.isGap.put(mode, new ArrayList<BitSet>());
                		}
                		
                		ArrayList<ArrayList<Integer>> t = this.recs.get(mode);
                		t.addAll(r.recs.get(mode));
                		
                		ArrayList<BitSet> g = this.isGap.get(mode);
                		g.addAll(r.isGap.get(mode));
                		
                	}     
                }
                
               
                
        
        		boolean isConsistent(ArrayList<ArrayList<Integer>> r, ArrayList<BitSet> g, int depth){
        		
        			for(int i=0; i<this.recs.get(modes.get(depth)).size(); i++){
        				ArrayList<ArrayList<Integer>> t = new ArrayList<ArrayList<Integer>>();
        				t.addAll(r);
        				t.add(this.recs.get(modes.get(depth)).get(i));
        				
        				ArrayList<BitSet> h= new ArrayList<BitSet>();
        				h.addAll(g);
        				h.add(this.isGap.get(modes.get(depth)).get(i));
        				
        				boolean is;
        				if(depth < modes.size()-1){
        					is = isConsistent(t,h,depth+1);
        				}else{
        					is = ic(t, h);
        				}
        				
        				if(is){
        					return true;
        				}
        			}
        			return false;
        		}
        		
        		
                private boolean ic(ArrayList<ArrayList<Integer>> recs, ArrayList<BitSet> isGap){
                	ArrayList<BitSet> recsinbitset = new ArrayList<BitSet>();
                	ArrayList<BitSet> isGapmask = new ArrayList<BitSet>();
                	
                	for(int i=0; i<recs.size(); i++){
                		ArrayList<Integer> rec = recs.get(i);
                		BitSet r = new BitSet();
                		BitSet g = new BitSet();
                		
                		BitSet t = isGap.get(i);
                		
                		for(int j=0; j<rec.size(); j++){
                			int m = rec.get(j);
                			r.set(m);
                			if(t.get(j)){
                				if(j>0) 
	                				for(int k= rec.get(j-1) + 1; k<m;k++)
	                					g.set(k);
                				else 
                					for(int k= 1; k<m;k++)
	                					g.set(k);
                			}
                		}
                		recsinbitset.add(r);
                		isGapmask.add(g);
                	}
                	
                	
                	
                	for(int i=0; i<recsinbitset.get(0).length(); i++){               
                		Boolean ic = null;
                		for(int j=0; j<recsinbitset.size(); j++){
                			boolean c = recsinbitset.get(j).get(i);
                			if(isGapmask.get(j).get(i)) continue;
                			
                			if(ic == null) ic = c;
                			else if(ic.equals(c)) continue;
                			else {
                			//	System.out.println(recs);

                            //	System.out.println(isGap);
                            	
                            //	System.out.println(recsinbitset);
                            //	System.out.println(isGapmask);
                            //	System.out.println();
                				return false;
                			}
                		}
                	}
                	
                	return true;
                	
                }
                
                public int hashCode(){
                        return title.hashCode();
                }
                
                public boolean equals(Object obj){
                        if(! (obj instanceof reconstructions)) return false;
                        reconstructions o = (reconstructions) obj;
                        
                        return this.title.equals(o.title);
                        
                }
                
                public String toString(){
                        return title + " " + modes;
                }
                
        }
        
        private static HashSet<ArrayList<Integer>> getPeptides(String fasta) throws IOException{
        	 HashSet<ArrayList<Integer>> out = new HashSet<ArrayList<Integer>>();
        	 if(fasta.isEmpty()) return out;
        	 BufferedLineReader in = new BufferedLineReader(fasta);
                String s;
               
                String t = "";
                while(true){
                        s=in.readLine();
                        if((s == null || s.startsWith(">")) && !t.isEmpty()){
                                //System.out.println(t);
	                        	for(int o =1; o<=1;o++){
	                                for(int i=o; i<t.length(); i++){
	                                        if(t.charAt(i-o) == 'K' || t.charAt(i-o) == 'R'){//12910
	                                                String p = "";
	                                                
	                                                int missedcleavage = 0;
	                                                for(int j=i;j<t.length(); j++){
	                                                        p += t.charAt(j);
	                                                        if(p.length() > 30 || missedcleavage >= 2) break;
	                                                        if(p.length() < 5) continue;
	                                                        if(p.endsWith("K") || p.endsWith("R")){
	                                                        		missedcleavage ++;
	                                                                ArrayList<Integer> pep = new ArrayList<Integer>();
	                                                                int pre = 0;
	                                                        
	                                                                for(int k=0; k<p.length(); k++){
	                                                                        AminoAcid aa = aaSet.getAminoAcid(p.charAt(k));
	                                                                        if(aa == null) break;
	                                                                        pre += aa.getNominalMass();
	                                                                        pep.add(pre);
	                                                                }
	                                                                if(pep.size() == p.length()){
	                                                                        out.add(pep);
	                                                                  //    System.out.println(p + "\t" + pep);
	                                                                }
	                                                        }
	                                                        
	                                                }
	                                        }else continue;
	                                }
	                        	}
                                t = "";
                                continue;
                        }
                        if(s == null) break;
                        if(s.startsWith(">")) continue;
                        t += s;
                   //     if(out.size() > 500000)break;
                }
                System.out.println(fasta + "\t" + out.size());
                return out;
        }
        
        
        public static void main(String[] args) throws IOException {
                String[] files = {                             
                             "/home/kwj/Dropbox/HCDETDnew_1_2_5_7.txt",
                            // "/home/kwj/Dropbox/HCDETDnew_1_2_5_7.txt",
                           //  "/home/kwj/Dropbox/HCDETDnew_2_2_5_7.txt",
                             //"/home/kwj/Dropbox/HCDETD_2_2_20.txt",
                            // "/home/kwj/Dropbox/Zubarev_HCD_Annotated.mgf"
                };
  
                // target = "/home/kwj/workspace/inputs/DataBases/ipi.HUMAN.v3.52.decoy.fasta";

                //HashSet<ArrayList<Integer>> decoys = getPeptides(decoy);
           
                // float[] Noffset = new float[251];
                // float[] Coffset = new float[251];
                float[] Ioffset = new float[251];
                
                HashSet<String> peptides = new HashSet<String>();
                HashSet<String> peptidesWithPTM = new HashSet<String>();
                
                ArrayList<reconstructions> output2 = new ArrayList<reconstructions>();
                HashSet<reconstructions> output = new HashSet<reconstructions>();
                HashSet<String> titles = new HashSet<String>();
                
                for(int mode = files.length-1; mode >=0; mode --){
                        BufferedLineReader in = new BufferedLineReader(files[mode]);
                        
                        String s;
                        
                        String title=""; 
                //      int sn=0; 
                        ArrayList<String> recs =null;
                       
                        boolean matchmode = false;
                        boolean matched = false;
                        boolean PTMmatched = false;
                        
                        while((s=in.readLine())!= null){
                                if(mode == files.length-1 && files[mode].endsWith("mgf")){
                                        if(s.startsWith("TITLE")) title = s.substring(6);
                                        if(s.startsWith("SEQ")){
                                                recs = new ArrayList<String>();
                                                recs.add(s.substring(4));
                                        }
                                      //  if(s.startsWith("CHARGE=2+")){
                                        	titles.add(title);
                                                reconstructions r = new reconstructions(title, mode, recs);
                                               // r.hit = true;
                                                if(output.contains(r)){
                                                        output2.get(output2.indexOf(r)).add(r);
                                                }else{
                                                        output2.add(r);
                                                        output.add(r);
                                                }
                                      //  }
                                        continue;
                                }
                                
                                if(s.isEmpty()) break;
                                if(s.startsWith("Set Accuracy")){
                                        continue;
                                }
                                String[] token = s.split("\t");
                                
                                if(s.startsWith(">>")){
                                	
                                	if(recs != null && !recs.isEmpty()){
                                		reconstructions r = new reconstructions(title, mode, recs);
                                	    r.hit = matched; 
	                                    r.PTMhit = PTMmatched;
	                                    if(output.contains(r)){
	                                    	output2.get(output2.indexOf(r)).add(r);
	                                    }else{
	                                    	output2.add(r);
	                                    	output.add(r);
	                                    }
                                    }
                                	
                                    title = token[3];
                                  
                            //      sn = Integer.parseInt(token[2]);
                                    matchmode = false;
                                    matched = false;
                                    PTMmatched = false;
                                    recs = new ArrayList<String>();
                                }else if(s.startsWith("Matches:")){
                                	matchmode = true;
                                }else if(!matchmode){
                                    recs.add(token[0]);
                                }else{
                                	if(token.length == 1){
                                		
                                		if(!matched)
                                		{
                                			if(token[0].endsWith("K") || token[0].endsWith("R")){
                                			//	System.out.println(s+" "+title);
                                				if(!titles.contains(title)) 
                                					peptides.add(token[0]);  
                                				matched = true;
                                			}
                                		}                     		
                                	}
                                	else{
                                	//	float pm = new Peptide(token[0], aaSet).getMass() + Float.parseFloat(token[1]);
                            			float of = Float.parseFloat(token[1])+50;
                            			
                            		//	if(Float.parseFloat(token[2]) == 0 || Math.abs(Float.parseFloat(token[3]) - pm) < 1)
                                			
	                                		if(!PTMmatched){
	                                			if(token[0].endsWith("K") || token[0].endsWith("R")){
	                                				if(!titles.contains(title))
	                                					peptidesWithPTM.add(token[0]);
	                                				PTMmatched = true; 
	                                			}
	                                		}
                                		
                                		
                                		if(token[0].endsWith("K") || token[0].endsWith("R")){
                                			
                                			//if(Float.parseFloat(token[2]) == 0 || Math.abs(Float.parseFloat(token[3]) - pm) < 1)
                                			
                                				Ioffset[Math.round(of)]++;
                                			/*
                                			if(Float.parseFloat(token[2]) == 0){//N
                                				Noffset[Math.round(of)]++;
                                			}else if(Math.abs(Float.parseFloat(token[3]) - pm) < 1){//C
                                				Coffset[Math.round(of)]++;
                                				//System.out.println(recs + " " + token[0] + " " + (of-50));
                                			}else{//I
                                				Ioffset[Math.round(of)]++;
                                			} */                              			
                                			
                                		}
                                		
                                		
                                		
                                	}
                                }
                                
                        }
                        in.close();
                }
                
                HashMap<ArrayList<Integer>, int[]> result = new HashMap<ArrayList<Integer>, int[]>();
                
                ArrayList<Integer> p = new ArrayList<Integer>(); p.add(0);// p.add(1);
                p.add(2);

             //   int numpep = 0;
                //    int cn = 0;
                
                HashMap<Integer, ArrayList<reconstructions>> outputpm = new HashMap<Integer, ArrayList<reconstructions>>();
                
                for(reconstructions r : output){
                        
                	for(int key : r.recs.keySet()){
                		int y = r.recs.get(key).get(0).get(r.recs.get(key).get(0).size()-1);
                        if(!outputpm.containsKey(y))
                            outputpm.put(y, new ArrayList<reconstructions>());
                        
                        ArrayList<reconstructions> v = outputpm.get(y);
                        v.add(r);
                	}             
                }
                
                   /* for(ArrayList<Integer> pep : targets){
                         //   System.out.println(cn++);
                            boolean targethit = false;
                            
                            if(outputpm.containsKey(pep.get(pep.size()-1))){
                                    for(reconstructions r : outputpm.get(pep.get(pep.size()-1))){
                                        if(r.modes.contains(3)) continue;//TODO
                                    	for(int key : r.recs.keySet()){
                                            for(ArrayList<Integer> rec : r.recs.get(key)){
                                                    
                                                //if(!rec.get(rec.size()-1).equals(pep.get(pep.size()-1))) break;
                                                
                                                if(pep.containsAll(rec)){
                                                  //      System.out.println(rec + "\t" + pep + "\t" + r);
                                                        r.hit.put(key, true);
                                                        targethit = true;
                                                        break;
                                                }
                                            }
                                    	}
                                          //  if(targethit) break;
                                    }
                            }
                            if(targethit) numpep++;
                            
                    }*/
                    
                
                
                for(reconstructions r : output){
                //	if(r.modes.contains(1)) continue;
                	
                		
            		//if(r.hit) continue;
            		//if(r.PTMhit) continue;
            		
            		if(!result.containsKey(r.modes))
                           result.put(r.modes, new int[2]);

                    
                    if(r.modes.equals(p)){
                    //      System.out.println(r + "\t" + r.title);
                    }
                    
                    int[] k = result.get(r.modes);
                    
                    k[0] ++;
                    if(r.isConsistent(new ArrayList<ArrayList<Integer>>(), new ArrayList<BitSet>(), 0)){ 
                    	k[1] ++;
                    	
                    	
                    }else if(r.modes.contains(2) && r.modes.contains(1)){
                    	//System.out.println(r.recs);
                    	//System.out.println(r.isGap + "\n");
                    }
                 //   else System.out.println(r);
                }
                
                
                
             
                
                int cr = 0, cr2 = 0 ;
                for(reconstructions r : output){
                	if(r.modes.contains(1)) continue;
                	  if(r.hit){
                    	  cr ++;
                    	  System.out.println(r);
                	  }
                	  if(r.PTMhit) cr2 ++;
                        //else System.out.println(r);
                }
                        
                System.out.println("Ioffset = [");
                for(int i=0; i<Ioffset.length; i++){
                	System.out.println((i-50) + "\t" + Ioffset[i]);
                }
                System.out.println("];");
                
                /*
                System.out.println("Coffset = [");
                for(int i=0; i<Coffset.length; i++){
                	System.out.println((i-50) + "\t" + Coffset[i]);
                }
                System.out.println("];");
                
                System.out.println("Noffset = [");
                for(int i=0; i<Noffset.length; i++){
                	System.out.println((i-50) + "\t" + Noffset[i]);
                }
                System.out.println("];");
                */
                
                for(ArrayList<Integer> key : result.keySet()){
                        int[] v = result.get(key);
                        
                        System.out.println(key + "\t" + v[0] + "\t" + v[1] + "\t" + (float)(v[1])/v[0]*100);
                        
                }
             	 System.out.println(peptides);
                System.out.println(peptides.size() + "\t" + peptidesWithPTM.size() + "\t" + cr + "\t" + cr2 +  "\t" + output.size());
        }

}