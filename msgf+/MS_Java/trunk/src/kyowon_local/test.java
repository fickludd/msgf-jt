package kyowon_local;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import parser.BufferedLineReader;

public class test {
	 public static void main(String args[]) throws IOException {
		 	String activation = "CID";//ETD
		 	String enzyme = "Tryp";//LysN
		    for(int charge = 2;charge<5;charge++){
		    String filename = "/home/kwj/Desktop/" + activation + enzyme + "Charge"+charge+".txt";
		    String title = activation+"-"+enzyme + " charge " + charge;
		    String outfilename = filename + ".m";
		    BufferedLineReader in = new BufferedLineReader(filename);
		    BufferedWriter out = new BufferedWriter(new FileWriter(outfilename));
		    String s;
		    boolean read = false;
		    int numPartition = Integer.MAX_VALUE;
		    int n = 1;
		    ArrayList<String> ionNames = new ArrayList<String>();
		    out.write("close;\nclear;\n");
		    while((s=in.readLine())!=null){
		    	
		    	if(s.startsWith("#RankDistributions")){
		    		numPartition = Integer.parseInt(s.split("\t")[1]);
		    		read = true;
		    		n=1;
		    	}
	  	
		    	if(s.startsWith("#ErrorDistributions")) read = false;
		    	
		    	if(read){
		    		if(numPartition == 0){
		    			String[] token = s.split("\t");
		    			if(token[0].startsWith("noise")) continue;
		    			ionNames.add(token[0]);
		    			out.write("ionDist(:," +n + ")=[");
		    			System.out.print("ionDist(:," +(n++) + ")=[");
		    			for(int i=1; i< token.length-2;i++){
		    				out.write(token[i] + " ");
		    				System.out.print(token[i] + " ");
		    			}
		    			
		    			out.write("];\n");
		    			System.out.println("];");
		    		}


		    		if(s.startsWith("Partition")) numPartition--;
		    		
		    	}
		    	
		    }
		    ionNames.add("Unexplained");
		    out.write("\nionDist(:," + n + ")= 1-sum(ionDist');\n");
		    out.write("figure1 = figure('XVisual',...\n'0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');\n");
		    out.write("axes1 = axes('Parent',figure1,'FontSize',20);\n");
		    out.write("box('on');hold('all');\n");
		    out.write("plot1 = plot(ionDist);\n");
		    
		    for(int i=0;i<ionNames.size();i++){
		    	if(i<ionNames.size()/2)
		    		out.write("set(plot1("+(i+1)+"),'DisplayName','"+ionNames.get(i).replace('_', ',')+"','LineWidth',1);\n");
		    	else if(i<ionNames.size()-1)
		    		out.write("set(plot1("+(i+1)+"),'DisplayName','"+ionNames.get(i).replace('_', ',')+"','Marker','*','LineStyle','none');\n");
		    	else
		    		out.write("set(plot1("+(i+1)+"),'DisplayName','"+ionNames.get(i).replace('_', ',')+"','LineWidth',1,'Color',[0 0 0]);\n");
		    }
		    out.write("xlabel('Rank','FontSize',20);\n");
		    out.write("ylabel('Probability','FontSize',20);\n");
		    out.write("title({'"+title+"'},'FontSize',20,'FontName','helvetica');\n");
		    out.write("legend(axes1,'show');\n");
		    out.write("xlim([1 100]);\nylim([0 1])\n");
		    out.close();
		    in.close();
		    }
		  }
}
