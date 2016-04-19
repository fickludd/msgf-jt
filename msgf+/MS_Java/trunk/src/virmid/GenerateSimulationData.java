package virmid;

import java.util.Random;

public class GenerateSimulationData {

    private int[] generateNTGenoType(double[] accG) {
        int[] gs = new int[2];
        Random r = new Random();
        double rn = r.nextDouble();
        int j = 0;
        for (int i = 0; i < accG.length; i++) {
            if (accG[i] >= rn) {
                j = i;
                break;
            }
        }

        gs[0] = j / 3;
        gs[1] = j - 3 * gs[0];

        return gs;
    }

    private boolean generatePointGenoType(double alpha, int[] gs, double e) {
        boolean genoType = true;

        Random r = new Random();

        int gN = gs[0];
        int gT = gs[1];

        //System.out.println(j + "\t" + gN + "\t" + gT);

        int g = gT;
        //	boolean isNormal = false;

        if (r.nextDouble() < alpha) {
            g = gN;
            //	isNormal = true;
        }


        if (g == 1) {
            genoType = r.nextBoolean();
        } else if (g == 2) {
            genoType = false;
        }

        if (genoType) {
            if (r.nextDouble() < e) {
                genoType = false;
            }
        } else {
            if (r.nextDouble() < e / 3) {
                genoType = true;
            }
        }
        //System.out.print((genoType? "T" : "F")  + " " + (isNormal? "N" : "M") + " : " );

        return genoType;

    }

    private double[] generateAccumulatedG(double[] G) {
        double[] accG = new double[9];

        for (int i = 0; i < 8; i++) {
            accG[i] = G[i] + (i > 0 ? accG[i - 1] : 0);
        }
        accG[8] = 1;

        return accG;
    }

    void generateData(boolean[][] X, double[][] e, double alpha, double[] G, int sampleNumber) {

        Random r = new Random();

        double[] accG = generateAccumulatedG(G);

        for (int i = 0; i < X.length; i++) {
            X[i] = new boolean[40 + r.nextInt(20)];
            e[i] = new double[X[i].length];

            int gs[] = generateNTGenoType(accG);

            for (int j = 0; j < X[i].length; j++) {
                e[i][j] = Math.max(0.00001, r.nextGaussian() / 500);
                e[i][j] = Math.min(e[i][j], 0.99999);

                X[i][j] = generatePointGenoType(alpha, gs, e[i][j]); // should be fixed accross j

            }
            //	System.out.println();
        }
    }

    public static void main(String[] args) {
        GenerateSimulationData test = new GenerateSimulationData();
        double[] G = {
            0.7, 0.03, 0.025,
            0.003, 0.1, 0.02,
            0.025, 0.02, 0.05,};

        double alpha = 0.1;
        int sampleNumber = 100;

        boolean[][] X = new boolean[sampleNumber][];
        double[][] e = new double[sampleNumber][];

        test.generateData(X, e, alpha, G, sampleNumber);

        /*
         * for(int i=0; i<X.length; i++){ for(int j=0; j<X[i].length; j++){
         * System.out.print(X[i][j]+" "); } System.out.println(); }
         */
    }
}
