package tdaAnalysis;
public class FisherExact {
    private static final boolean DEBUG = false;
    private double[] f;
    int maxSize;


    /**
     * constructor for FisherExact table
     *
     * @param maxSize is the maximum sum that will be encountered by the table (a+b+c+d)
     */
    public FisherExact(int maxSize) {
        this.maxSize = maxSize;
        double cf = 1.0;
        f = new double[maxSize + 1];
        f[0] = 0.0;
        for (int i = 1; i <= this.maxSize; i++) {
            f[i] = f[i - 1] + Math.log(i);
        }
    }

    /**
     * calculates the P-value for this specific state
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return the P-value
     */
    public final double getP(int a, int b, int c, int d) {
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p;
        p = (f[a + b] + f[c + d] + f[a + c] + f[b + d]) - (f[a] + f[b] + f[c] + f[d] + f[n]);
        return Math.exp(p);
    }

    /**
     * Calculates the one-tail P-value for the Fisher Exact test.  Determines whether to calculate the right- or left-
     * tail, thereby always returning the smallest p-value.
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return one-tailed P-value (right or left, whichever is smallest)
     */
    public final double getCumlativeP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;

        p += getP(a, b, c, d);
        if (DEBUG) {System.out.println("p = " + p);}
        if ((a * d) >= (b * c)) {
            if (DEBUG) {System.out.println("doing R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
            min = (c < b) ? c : b;
            for (i = 0; i < min; i++) {
                if (DEBUG) {System.out.print("doing round " + i);}
                p += getP(++a, --b, --c, ++d);
                if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
            }
            System.out.println("");
        }
        if ((a * d) < (b * c)) {
            if (DEBUG) {System.out.println("doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
            min = (a < d) ? a : d;
            for (i = 0; i < min; i++) {
                if (DEBUG) {System.out.print("doing round " + i);}
                double pTemp = getP(--a, ++b, ++c, --d);
                if (DEBUG) {System.out.print("\tpTemp = " + pTemp);}
                p += pTemp;
                if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
            }
        }
        return p;
    }

    /**
     * Calculates the right-tail P-value for the Fisher Exact test.
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return one-tailed P-value (right-tail)
     */
    public final double getRightTailedP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;

        p += getP(a, b, c, d);
        if (DEBUG) {System.out.println("p = " + p);}
        if (DEBUG) {System.out.println("doing R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
        min = (c < b) ? c : b;
        for (i = 0; i < min; i++) {
            p += getP(++a, --b, --c, ++d);

        }
        return p;
    }

    /**
     * Calculates the left-tail P-value for the Fisher Exact test.
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return one-tailed P-value (left-tail)
     */
    public final double getLeftTailedP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;

        p += getP(a, b, c, d);
        if (DEBUG) {System.out.println("p = " + p);}
        if (DEBUG) {System.out.println("doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
        min = (a < d) ? a : d;
        for (i = 0; i < min; i++) {
            if (DEBUG) {System.out.print("doing round " + i);}
            double pTemp = getP(--a, ++b, ++c, --d);
            if (DEBUG) {System.out.print("\tpTemp = " + pTemp);}
            p += pTemp;
            if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
        }


        return p;
    }


    /**
     *   Calculates the two-tailed P-value for the Fisher Exact test.
     *
     *   In order for a table under consideration to have its p-value included
     *   in the final result, it must have a p-value less than the original table's P-value, i.e.
     *   Fisher's exact test computes the probability, given the observed marginal
     *   frequencies, of obtaining exactly the frequencies observed and any configuration more extreme.
     *   By "more extreme," we mean any configuration (given observed marginals) with a smaller probability of
     *   occurrence in the same direction (one-tailed) or in both directions (two-tailed).
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return two-tailed P-value
     */
    public final double getTwoTailedP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;

        double baseP = getP(a, b, c, d);
//         in order for a table under consideration to have its p-value included
//         in the final result, it must have a p-value less than the baseP, i.e.
//         Fisher's exact test computes the probability, given the observed marginal
//         frequencies, of obtaining exactly the frequencies observed and any configuration more extreme.
//         By "more extreme," we mean any configuration (given observed marginals) with a smaller probability of
//         occurrence in the same direction (one-tailed) or in both directions (two-tailed).

        if (DEBUG) {System.out.println("baseP = " + baseP);}
        int initialA = a, initialB = b, initialC = c, initialD = d;
        p += baseP;
        if (DEBUG) {System.out.println("p = " + p);}
        if (DEBUG) {System.out.println("Starting with R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
        min = (c < b) ? c : b;
        for (i = 0; i < min; i++) {
            if (DEBUG) {System.out.print("doing round " + i);}
            double tempP = getP(++a, --b, --c, ++d);
            if (tempP <= baseP) {
                if (DEBUG) {System.out.print("\ttempP (" + tempP + ") is less than baseP (" + baseP + ")");}
                p += tempP;
            }
            if (DEBUG) {System.out.println(" a=" + a + " b=" + b + " c=" + c + " d=" + d);}
        }

        // reset the values to their original so we can repeat this process for the other side
        a = initialA;
        b = initialB;
        c = initialC;
        d = initialD;

        if (DEBUG) {System.out.println("Now doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
        min = (a < d) ? a : d;
        if (DEBUG) {System.out.println("min = " + min);}
        for (i = 0; i < min; i++) {
            if (DEBUG) {System.out.print("doing round " + i);}
            double pTemp = getP(--a, ++b, ++c, --d);
            if (DEBUG) {System.out.println("  pTemp = " + pTemp);}
            if (pTemp <= baseP) {
                if (DEBUG) {System.out.print("\ttempP (" + pTemp + ") is less than baseP (" + baseP + ")");}
                p += pTemp;
            }
            if (DEBUG) {System.out.println(" a=" + a + " b=" + b + " c=" + c + " d=" + d);}
        }
        return p;
    }
//
//    public static void main(String[] args) {
//
//        if(args.length != 4){
//            System.out.println("Please enter 4 values");
//            System.exit(0);
//        }
//        int[] argInts = new int[args.length];
//
//        for(int i = 0; i < argInts.length; i++){
//            argInts[i] = Integer.parseInt(args[i]);
//        }
//        FisherExact fe = new FisherExact(100);
//
//        System.out.println("\n*****Original algorithm");
//        double cumulativeP = fe.getCumlativeP(argInts[0], argInts[1], argInts[2], argInts[3]);
//        System.out.println("cumulativeP = " + cumulativeP );
//
//        System.out.println("\n*****Modified algorithm");
//        double algorithmSelectedP = fe.getAlgorithmSelected(argInts[0], argInts[1], argInts[2], argInts[3]);
//        System.out.println("algorithmSelectedP = " + algorithmSelectedP);
//
//        System.out.println("\n*****Left Tailed");
//        double leftTailedP = fe.getLeftTailedP(argInts[0], argInts[1], argInts[2], argInts[3]);
//        System.out.println("leftTailedP = " + leftTailedP);
//
//        System.out.println("\n*****Right Tailed");
//        double rightTailedP = fe.getRightTailedP(argInts[0], argInts[1], argInts[2], argInts[3]);
//        System.out.println("rightTailedP = " + rightTailedP);
//
//        System.out.println("\n*****Two Tailed");
//        double twoTailedP = fe.getTwoTailedP(argInts[0], argInts[1], argInts[2], argInts[3]);
//        System.out.println("twoTailedP = " + twoTailedP);
//    }

    
    static int[] getTable(int t1, double f1, int t2, double f2){
    	int[] ret ={
    			t1,
    			(int) Math.round(t1*f1),
    			t2,
    			(int) Math.round(t2*f2),	
    	};
    	return ret;
    	
    }
    
    public static void main(String[] args) {

        int[][] argInts = new int[10][4];
       
       argInts[0] =  getTable(1835,0.03,1909,0.037);
       argInts[1] =  getTable(1157,0.013,1136,0.010);
       argInts[2] =  getTable(986,0.005,1078,0.013);
       argInts[3] =  getTable(781,0.008,869,0.016);
       
       argInts[4] =  getTable(1835,0.03,1749,0.021);
       argInts[5] =  getTable(986,0.005,902,0.004);
       argInts[6] =  getTable(781,0.008,781,0.008);
       
       argInts[7] =  getTable(1835,0.03,1786,0.029);
       argInts[8] =  getTable(1157,0.013,986,0.005);
       argInts[9] =  getTable(1157,0.013,911,0.014);
       
      // argInts[10] =  getTable(781,0.008,781,0.008);
      // argInts[11] =  getTable(781,0.008,781,0.008);
       
        FisherExact fe = new FisherExact(10000);

        for (int i = 0; i < argInts.length; i++) {
            System.out.println("\na=" + argInts[i][0] + " b=" + argInts[i][1] + " c=" + argInts[i][2] + " d=" + argInts[i][3]);
            System.out.print("*****Original algorithm: ");
            double cumulativeP = fe.getCumlativeP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("\tcumulativeP = " + String.format("%.1f %%", cumulativeP*100));

           // System.out.print("*****Two Tailed: ");
           // double twoTailedP = fe.getTwoTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
          //  System.out.println("\ttwoTailedP = " + twoTailedP);
            /*System.out.print("*****Left Tailed: ");
            double leftTailedP = fe.getLeftTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("\tleftTailedP = " + leftTailedP);

            System.out.print("*****Right Tailed: ");
            double rightTailedP = fe.getRightTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("\trightTailedP = " + rightTailedP);

            System.out.print("*****Two Tailed: ");
            double twoTailedP = fe.getTwoTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("\ttwoTailedP = " + twoTailedP);*/
        }
    }
}