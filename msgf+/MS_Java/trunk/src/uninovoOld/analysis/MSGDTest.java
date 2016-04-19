package uninovoOld.analysis;

import java.io.IOException;
import java.util.ArrayList;

import msgf.NominalMass;
import msutil.Peptide;
import parser.BufferedLineReader;

public class MSGDTest {
	static public void main(String[] args) throws IOException{
		String file = "/home/kwj/Dropbox/CID_Tryp_Confident_test_MSGD.grc";
		int charge = 2;
		
		int sn=0, csn = 0, c = 0;
		
		BufferedLineReader in = new BufferedLineReader(file);
		String s;
		ArrayList<Integer> correctNominalMasses = null;
		ArrayList<Integer> gapMasses = null;
		boolean isCorrect = true;
		Peptide pep = null;
		
		while((s=in.readLine())!=null){
			if(s.startsWith("#")){
				if((c = Integer.parseInt(s.split("\t")[4])) != charge) continue;
				
				sn ++;
				
				if(!isCorrect){
					System.out.println(pep);
				}
				
				correctNominalMasses = new ArrayList<Integer>();
				pep = new Peptide(s.split("\t")[6]);
				
				for(float m : pep.getPRMMasses(true, 0)){
					correctNominalMasses.add(NominalMass.toNominalMass(m));
				}
				
				correctNominalMasses.add(pep.getNominalMass());
				isCorrect = false;
				
				continue;
			}
			if(isCorrect || c != charge) continue;
			
			gapMasses = new ArrayList<Integer>();
			String[] t = s.substring(1, s.length()-1).split(",");
			int pm = 0;
			
			for(String u : t){
				pm += Integer.parseInt(u.trim());
				gapMasses.add(pm);
			}
			if(correctNominalMasses.containsAll(gapMasses)){
				csn++;
				//System.out.println(correctNominalMasses + "\t" + gapMasses);
				isCorrect = true;
			}
			
			
		}
		
		System.out.println(sn + "\t" + csn);
		
		
		in.close();
	}
}
