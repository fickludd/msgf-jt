package hybridsearch;

import java.io.IOException;
import java.io.PrintStream;

import parser.BufferedLineReader;

public class AnnotationChanger {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		BufferedLineReader in = new BufferedLineReader(args[0]);
		PrintStream out = new PrintStream(args[1]);
		boolean isDecoy = args[2].equals("1"); 
		
		String s;
		
		while((s=in.readLine())!=null){
			String[] t = s.split("\t");
			
			if(s.startsWith("Scan#")){
				for(int a=0;a<16;a++){//TODO erase..
					out.print(t[a]+"\t");
				}
				out.print(t[16]);
				out.println();
				continue;
			}
			
			
			//if(t[6].contains("C")){
			t[2] = t[2].replace("C", "C-57.02146");
			//}
			
			if(t[2].contains("(")){
				String tt = new String(t[2]);
				while(true){
					int j = tt.indexOf("(");
					int i = tt.indexOf(")");
					if(i>0){
						String[] k = tt.substring(j+1, i).split(",");
						String aa = k[0];
						float off = Float.parseFloat(k[1]);
						if(aa.equals("C-57.02146")){
							off -= 57.02146f;
							aa = "C";
						}
						
						String rep = aa;
						if(off != 0){
							rep += (off>0? "+" : "") + off;
						}
						tt = tt.replace(tt.substring(j, i+1), rep);
					
					}else break;					
				}
				
				t[2] = tt;
				
				
			}
			
			//if((Math.abs(new Peptide(t[2]).getParentMass() -18) -	(Float.parseFloat(t[2])-1) * Integer.parseInt(t[9]))>1  ){
			//	System.out.println(t[2] + "\t" + s.split("\t")[2]);
			//}
			if(isDecoy) t[12] = "1";
			
			for(int a=0;a<t.length;a++){
				out.print(t[a]+"\t");
			}
			out.println();
			
		}
		
		out.close();
		in.close();
	}

}
