package uninovoOld;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import parser.BufferedLineReader;

public class ScoringDataIterator  implements  Iterator<ScoringData>{
	private BufferedLineReader reader;
	private AminoAcidSet aaSet;
	private String nextChunk = "";
	private Enzyme enzyme;
	private ArrayList<Tolerance> tols = null;
	private Tolerance pmtol = null;
	
	public ScoringDataIterator(String file, AminoAcidSet aaSet, Enzyme enzyme){
		this.aaSet = aaSet;
		this.enzyme = enzyme;
		try {
			reader = new BufferedLineReader(file);
			getNextChunk();
		} catch (IOException ioe) {
			 System.err.println(ioe);
		     System.exit(-9);
		}
	}
	
	private void getNextChunk(){
		if(nextChunk.endsWith("*#*\n")){
			nextChunk = "";
			return;
		}
		nextChunk = "";
		String s;
		while((s=reader.readLine())!=null){
			if(s.startsWith(">>TOL")){
				if(tols == null) tols = new ArrayList<Tolerance>();
				tols.add(Tolerance.parseToleranceStr(s.split("\t")[1]));
				continue;
			}else if(s.startsWith(">>PMTOL")){
				pmtol = Tolerance.parseToleranceStr(s.split("\t")[1]);
				continue;
			}
			
			nextChunk+=s+"\n";
			
			if(s.startsWith("***")){
				break;
			}else if(s.startsWith("*#*")){
				break;
			}
			
			
			
		}
	}
	
	
	@Override
	public boolean hasNext() {
		return !nextChunk.isEmpty();
	}

	@Override
	public ScoringData next() {
		if(hasNext()){
			int specID = 0;
			int[] charges = null;
			float pm = 0;
			String specfilename = "";
			
			ArrayList<Node> nodes = new ArrayList<Node>();
			SpectrumGraph graph = null;
			ArrayList<DenovoReconstruction> recs = new ArrayList<DenovoReconstruction>();
			
			for(String s : nextChunk.split("\n")){
				if(s==null || s.isEmpty() || s.startsWith("*")) continue;
				if(s.startsWith("#")){
					String[] token = s.split("\t");
					specID = Integer.parseInt(token[0].substring(1));
					String[] cs = token[4].split(",");
					charges = new int[cs.length];
					for(int l=0; l<charges.length;l++){
						charges[l] = Integer.parseInt(cs[l]);
					}

					pm = Float.parseFloat(token[3]);
					specfilename = token[2];
				}else if(s.startsWith(":")){
					String[] token = s.substring(1, s.length()-1).split("::");

					for(String t : token){						
						String[] tk = t.substring(0, t.length()).split(",");
						float mass = Float.parseFloat(tk[0]);
						float score = Float.parseFloat(tk[2]);
						Node node = new Node(mass, tols.get(Integer.parseInt(tk[1])), pm, Integer.parseInt(tk[1]));
						node.setLR(score);
						nodes.add(node);
					}		
					graph = new SpectrumGraph(nodes, enzyme, pmtol, aaSet, charges, tols.size(), false);
				}else{
					ArrayList<Integer> indices = new ArrayList<Integer>();
					String[] token = s.split(",");
					for(String t : token){
						indices.add(Integer.parseInt(t));
					}
					
					recs.add(new DenovoReconstruction(graph, indices));
				}
			}
			
			ScoringData sd = new ScoringData(specID, graph, recs, specfilename);
			getNextChunk();
			
			return sd;
		}else return null;
	}

	@Override
	public void remove() {
		return;		
	}

}
