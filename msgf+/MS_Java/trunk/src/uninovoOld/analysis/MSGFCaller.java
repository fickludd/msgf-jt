package uninovoOld.analysis;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

public class MSGFCaller {
	
     public static void main(String args[]) throws IOException {

       if (args.length <= 0) {
         System.err.println("Need command to run");
         System.exit(-1);
       }

       Process process = new ProcessBuilder(args).start();
       InputStream is = process.getInputStream();
       InputStreamReader isr = new InputStreamReader(is);
       BufferedReader br = new BufferedReader(isr);
       String line;

       System.out.printf("Output of running %s is:", 
          Arrays.toString(args));

       while ((line = br.readLine()) != null) {
         System.out.println(line);
       }

     }
		 
}
