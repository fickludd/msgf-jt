package virmid;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Random;

import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.linear.LinearConstraint;
import org.apache.commons.math3.optimization.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optimization.linear.Relationship;
import org.apache.commons.math3.optimization.linear.SimplexSolver;
import org.apache.commons.math3.util.ArithmeticUtils;

public class MLE {

    private boolean estAlpha;
    private boolean plotForMatlab = false;
    private double alpha = -1;
    private double[][][][] f = null;
    private double[][][] af = null;
    private double[][] logx = null;
    // private int[] depth = null;
    //private float[] allAs = null;
    private double[] V = null;
    private int sc = 0;
    private double x = 0;
    private double y = 0;
    private double beta = 0.99f;
    //private int nSample = 0;
    private double delta = 1e-100;
    private double minBratio = -1;
    private double maxBratio = -1;// 0.25;
    // private HashMap<ArrayList<Integer>, Float> coefs;

    // nSample should exceed the actual sample number!
    public MLE(int nSample, double alpha) { // genotype est or calling
        this.alpha = alpha;
        this.estAlpha = false;
        this.logx = new double[nSample][9];
    }

    public MLE(int nSample, double minBratio, double maxBratio) { // alpha estimation
        this.estAlpha = true;
        this.minBratio = minBratio;
        this.maxBratio = maxBratio;

        // double t1 = (1 - y) * (1 + alpha);
        // double t2 = (1 - x) * (1 - alpha);
        //this.alpha = (t1 - t2) / (t1 + t2);
        this.f = new double[2][nSample][][];
        this.af = new double[nSample][][];
    }


    /*
     * private float getLogDenominator(int k, int l, float alpha, int d, int i,
     * int n){ float p = 0; float deno = 0;
     *
     * if(k==0){
     *
     * }else if(k==1){ p += alpha * .5; }else if(k==2){ p += alpha; }
     *
     * if(l==0){
     *
     * }else if(l==1){ p += (1-alpha) * .5; }else if(l==2){ p += 1-alpha; }
     *
     * int start = 1; if(p == 0){ p = 0.01f; deno = allAs[i]; start = 1; }
     *
     * for(int j=start;j<n;j++){ ArrayList<Integer> key = new
     * ArrayList<Integer>(); key.add(d); key.add(j); Float coef =
     * coefs.get(key); if(coef == null){ coef =
     * (float)ArithmeticUtils.binomialCoefficient(d, j); coefs.put(key, coef); }
     * deno += coef * Math.pow(p, j) * Math.pow(1 - p, d - j); } deno = 1-deno;
     * deno = (float)Math.log(deno);
     *
     * return deno; }
     *
     *
     *
     * private float getLogDenominator(int k, int l, float alpha, int d, float
     * allA, int n){ float p = 0; float deno = 0;
     *
     * if(k==0){
     *
     * }else if(k==1){ p += alpha * .5; }else if(k==2){ p += alpha; }
     *
     * if(l==0){
     *
     * }else if(l==1){ p += (1-alpha) * .5; }else if(l==2){ p += 1-alpha; }
     *
     * int start = 1; if(p == 0){ p = 0.01f; deno = allA; start = 1; }
     *
     * for(int j=start;j<n;j++){ ArrayList<Integer> key = new
     * ArrayList<Integer>(); key.add(d); key.add(j); Float coef =
     * coefs.get(key); if(coef == null){ coef =
     * (float)ArithmeticUtils.binomialCoefficient(d, j); coefs.put(key, coef); }
     * deno += coef * Math.pow(p, j) * Math.pow(1 - p, d - j); } deno = 1-deno;
     * deno = (float)Math.log(deno);
     *
     * return deno; }
     */
    private double getPValue(double[][] af, double alpha, double minScore, int k, int l) {
        double pVal = 0;

        double[][] d = new double[af.length + 1][(int) minScore + 5];
        HashSet<Integer> nonZeros = new HashSet<Integer>();

        d[0][0] = 1;
        nonZeros.add(0);

        for (int i = 1; i < d.length; i++) {
            double ha = h(alpha, af[i - 1][k], af[i - 1][l], k, l);
            int sa = (int) (-Math.log(ha) + 0.5);
            int sb = (int) (-Math.log(1 - ha) + 0.5);
            HashSet<Integer> newNonZeros = new HashSet<Integer>();

            for (int nz : nonZeros) {
                int n1 = nz + sa;
                if (n1 >= 0 && n1 <= minScore) {
                    d[i][n1] += d[i - 1][nz] * ha;
                    newNonZeros.add(n1);
                }
                int n2 = nz + sb;
                if (n2 >= 0 && n2 <= minScore) {
                    d[i][n2] += d[i - 1][nz] * (1 - ha);
                    newNonZeros.add(n2);
                }
            }
            nonZeros = newNonZeros;
        }

        for (int i : nonZeros) {
            if (d[d.length - 1][i] == 0) {
                continue;
            }

            //System.out.println(d.length + "\t" + minScore + "\t" + i + "\t" + d[d.length - 1][i]);

            pVal += d[d.length - 1][i];// * Math.pow(Math.E, -i);
        }
        //System.out.println("***" + pVal);
        return pVal;
    }

    private double[] getDenoTable(double[][] af, double alpha, double minBratio, double maxBratio) {
        double[] denoTable = new double[9];

        int minBnum = (int) (Math.ceil(af.length * minBratio) + 0.01);//-(int) (-af.length * minBratio);
        minBnum = Math.max(1, minBnum);


        int maxBnum = (int) (Math.ceil(af.length * maxBratio) + 0.01);
        maxBnum = Math.max(minBnum + 1, maxBnum);
        //System.out.println(af.length + "\t" + minBnum);
        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                double[][] d = new double[af.length + 1][maxBnum];

                d[0][0] = 1e40;

                for (int i = 1; i < d.length; i++) {
                    double fk = af[i - 1][k];
                    double fl = af[i - 1][l];

                    double ha = h(alpha, fk, fl, k, l);

                    for (int j = 0; j < maxBnum; j++) {
                        if (d[i - 1][j] == 0 && (j == 0 || d[i - 1][j - 1] == 0)) {
                            continue;
                        }
                        d[i][j] = d[i - 1][j] * ha;
                        if (j > 0) {
                            double hb = 1 - ha;//h(alpha, 1 - af[i -1][k], 1 - af[i - 1][l], k, l);
                            d[i][j] += d[i - 1][j - 1] * hb;
                        }
                    }
                }

                double v = 0;//d[0][0];

                if (maxBratio > 0) {
                    for (int i = minBnum; i < maxBnum; i++) {
                        v += d[af.length][i];
                    }
                } else {
                    v = d[0][0];
                    for (int i = 0; i < minBnum; i++) {
                        v -= d[af.length][i];
                    }
                }
                // for (int i = 0; i < minBnum; i++) {
                //      v -= d[af.length][i];
                // }



                //System.out.print(v+" ");
                v = Math.max(v, d[0][0] * 1e-40f);

                denoTable[k * 3 + l] = (d[0][0] / v);
                //System.out.println(denoTable[k * 3 + l]+" ");
            }
        }
        // System.out.println();
        //MC.printVector(denoTable);
        return denoTable;
    }

    /*
     * private float[] getDenoTable(float[][] af, float[][] f, float alpha, int
     * minBnum, float minBratio) { float[] denoTable = new float[9]; if
     * (minBratio >= 0) { minBnum = (int) (af.length * minBratio + 1); minBnum =
     * Math.max(1, minBnum); // if(estAlpha) System.out.println(af.length + "\t"
     * + ((float)af.length * minBratio) + "\t" + minBratio + "\t" + minBnum); }
     *
     *
     * //System.out.println(af.length + "\t" + minBnum); for (int k = 0; k < 3;
     * k++) { for (int l = 0; l < 3; l++) {
     *
     * double[][] d = new double[af.length + 1][minBnum];
     *
     * d[0][0] = Double.MAX_VALUE;
     *
     * for (int i = 1; i < d.length; i++) { for (int j = 0; j < minBnum; j++) {
     * double ha = h(alpha, af[i - 1][k], af[i - 1][l], k, l);
     *
     * d[i][j] = d[i - 1][j] * ha; if (j > 0) { double hb = 1 - ha;//h(alpha, 1
     * - af[i - 1][k], 1 - af[i - 1][l], k, l); d[i][j] += d[i - 1][j - 1] * hb;
     * } } }
     *
     * double v = d[0][0];
     *
     * for (int i = 0; i < minBnum; i++) { v -= d[af.length][i]; }
     * //System.out.print(v+" "); //v = Math.max(v, 1e-10f);
     *
     * double s = 1; int diff = 0; for (int j = 0; j < f.length; j++) { float fk
     * = f[j][k]; float fl = f[j][l]; if (fk != af[j][k]) { diff++; } s *=
     * (h(alpha, fk, fl, k, l)); }
     *
     *
     *
     * if (v / d[0][0] < s) { System.out.println("SHIT " + minBratio + "\t" +
     * af.length + "\t" + diff + "\t" + minBnum + "\t" + v / d[0][0] + "\t" +
     * s); } denoTable[k * 3 + l] = (float) (d[0][0] / v);
     * //System.out.println(denoTable[k * 3 + l]+" "); } } //
     * System.out.println(); //MC.printVector(denoTable); return denoTable; }
     */
    private double f(int g, float er, float em, boolean isA, double beta) {
        //float x = 0.05f;
        double a;
        //float y = 0;
        double z = (1 - x) / (4 - 3 * y - x);//(4 - 3 * y - x);
        double a1 = z * er + (1 - er) * (1 - em + beta * em);//
        double a2 = z * er + beta * em * (1 - er);//((1-x) / (4 - 3 * y - x))
        if (g == 0) // aa
        {
            a = a1;
        } else if (g == 1) { // ab
            a = .5 * a1 + .5 * a2;
        } else { // bb
            a = a2;
        }
        a = isA ? a : 1 - a;

        return Math.max(a, Double.MIN_VALUE);
    }

    private double f(int g, float er, float em, boolean isA, double beta, double offset) {
        //float x = 0.05f;

        em = (float) Math.min(1, offset * em);

        /*
         * double p = f(g, er, em, true, beta);
         *
         * p-=offset; p = Math.max(p, Double.MIN_VALUE);
         *
         * if(!isA) p = 1-p; p = Math.max(p, Double.MIN_VALUE);
         */
        return f(g, er, em, isA, beta);
    }

    private double h(double alpha, double fk, double fl, int k, int l) {
        //float x = 0.05f;
        //alpha = 0.01f;
        double[] m = new double[4];
        if (k == 0) {
            m[0] = m[1] = alpha * (1 - x);
        } else if (k == 1) {
            m[0] = alpha * (1 - x);
            m[1] = alpha * (1 - y);
        } else {
            m[0] = m[1] = alpha * (1 - y);
        }

        if (l == 0) {
            m[2] = m[3] = (1 - alpha) * (1 - x);
        } else if (l == 1) {
            m[2] = (1 - alpha) * (1 - x);
            m[3] = (1 - alpha) * (1 - y);
        } else {
            m[2] = m[3] = (1 - alpha) * (1 - y);
        }

        double mul = (m[0] + m[1]) / (m[0] + m[1] + m[2] + m[3]);

        //if(k==0 && l == 1)
        //   System.out.println(mul + "\t" + alpha + "\t" + x + "\t" + m[0] + "\t"+ m[1] + "\t"+ m[2] + "\t"+ m[3] + "\t");

        double p = mul * fk + (1 - mul) * fl;

        // System.out.println(p + "\t" + (alpha * fk + (1 - alpha) * fl));
        return Math.max(p, Double.MIN_VALUE);
    }

    private double getLikelihoodAt(double[] V, int i) {
        double t = 0;
        if (estAlpha) {
            double[] denoTable = getDenoTable(af[i], V[0], minBratio, maxBratio);

            for (int k = 0; k < 3; k++) {
                if (V[k + 1] == 0) {
                    continue;
                }
                double u = 0;
                for (int j = 0; j < f[0][i].length; j++) {
                    u += Math.log(f[0][i][j][k]);
                }

                for (int l = 0; l < 3; l++) {
                    if (V[k * 3 + l + 4] == 0) {
                        continue;
                    }
                    double s = u;

                    for (int j = 0; j < f[1][i].length; j++) {
                        double fk = f[1][i][j][k];
                        double fl = f[1][i][j][l];
                        s += Math.log(h(V[0], fk, fl, k, l));
                    }

                    s += Math.log(denoTable[3 * k + l]);
                    t += Math.pow(Math.E, s) * V[k * 3 + l + 4] * V[k + 1];

                }
            }
        } else {
            for (int k = 0; k < 3; k++) {
                if (V[k + 1] == 0) {
                    continue;
                }
                for (int l = 0; l < 3; l++) {
                    if (V[k * 3 + l + 4] == 0) {
                        continue;
                    }
                    //int index = ;
                    //double s = (Math.log() + Math.log(V[k + 1]) + );

                    t += Math.pow(Math.E, logx[i][3 * k + l]) * V[k * 3 + l + 4] * V[k + 1];
                }
            }
        }

        if (t <= 0) {
            t = Double.MIN_VALUE;
        }
        return t;
    }

    private double[][] getCovarianceMatrix(double[] V) {
        double[][] c;
        double[][] G = new double[sc][3];

        for (int i = 0; i < sc; i++) {
            for (int v = 0; v < 5; v++) {
                double[] g = getGradientAt(V, i, v);
                if (v == 0) {
                    G[i][0] += g[0];
                }
                if (v == 2) {
                    G[i][1] += g[4];
                    G[i][2] += g[5];
                }
            }
        }

        // MC.printMatrix(MC.multiply(MC.transpose(G), G));
        c = MC.invert(MC.multiply(MC.transpose(G), G));

        return c;
    }

    private ArrayList<Integer> getBindingRaws(double[] V, int[][] A, int[][] b) {
        ArrayList<Integer> raws = new ArrayList<Integer>();

        for (int j = 0; j < A.length; j++) {
            double sum = 0;
            for (int i = 0; i < A[j].length; i++) {
                if (A[j][i] != 0) {
                    sum += A[j][i] * V[i];
                }
            }
            if (sum >= (double) b[j][0]) {
                raws.add(j);
            }
        }
        //   System.out.println(raws);
        return raws;
    }

    private double[][] getA1(int[][] A, ArrayList<Integer> br) {
        double[][] A1 = new double[br.size()][];

        int i = 0;
        for (int r : br) {

            A1[i] = new double[A[r].length];

            for (int j = 0; j < A1[i].length; j++) {
                A1[i][j] = A[r][j];
            }

            i++;

        }
        return A1;
    }

    private double[][] getAb2(int[][] Ab, ArrayList<Integer> br) {
        double[][] Ab2 = new double[Ab.length - br.size()][];

        int i = 0;
        for (int r = 0; r < Ab.length; r++) {
            if (br.contains(r)) {
                continue;
            }
            Ab2[i] = new double[Ab[r].length];
            for (int j = 0; j < Ab2[i].length; j++) {
                Ab2[i][j] = (double) Ab[r][j];
            }

            i++;
        }
        return Ab2;
    }

    private int[][] getA() {

        int[][][] As = {
            {
                {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// alpha >= 0
                {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// alpha <= 1
                //{0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// gNaa >= 0
                {0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,},// gTaa|gNaa >= 0
                {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,},// gTab|gNaa >= 0
                {0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,},// gTbb|gNaa >= 0
                // {0, 0, 0, 0, 10, -1, 0, 0, 0, 0, 0, 0, 0,},// gTab|gNaa >= gTaa|gNaa
                {0, 0, 0, 0, 0, -1, 1000, 0, 0, 0, 0, 0, 0,},// gTab|gNaa >= 10*gTbb|gNaa
            },
            {
                //    {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},//alpha >=0
                //   {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// alpha <=1
                {0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// gNaa >=0
                {0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// gNab >=0
                {0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// gNbb >=0
                {0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,},// gTaa|gNaa >=0
                {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,},// gTab|gNaa >=0
                {0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,},// gTbb|gNaa >=0
                {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0,},// gTaa|gNab >=0
                {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,},// gTab|gNab >=0
                {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0,},// gTbb|gNab >=0
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,},// gTaa|gNbb >=0
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0,},// gTab|gNbb >=0
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,},// gTbb|gNbb >=0
                {0, -1, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// gNaa >= 100*gNab 
                {0, 0, -1, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0,},// gNab >= 100*gNbb
                {0, 0, 0, 0, -1, 10000, 0, 0, 0, 0, 0, 0, 0,},// gTaa|gNaa >= 10000*gTab|gNaa
                //{0, 0, 0, 0, -1, 0, 10000, 0, 0, 0, 0, 0, 0,},// gTaa|gNaa >= 10000*gTbb|gTaa
                {0, 0, 0, 0, 0, -1, 100, 0, 0, 0, 0, 0, 0,},// gTab|gNaa >= 100*gTbb|gNaa
                {0, 0, 0, 0, 0, 0, 0, 1000000, -1, 0, 0, 0, 0,},// gTab|gNab >= 1000000*gTaa|gNab
                {0, 0, 0, 0, 0, 0, 0, 0, -1, 100, 0, 0, 0,},// gTab|gNab >= 100*gTbb|gNab
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000000, 0, -1,},// gTbb|gNbb >= 10000*gTaa|gNbb
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000000, -1,},// gTbb|gNbb >= 10000*gTab|gNbb
            //{0, 0, 0, 0, 0, 0, 0, 100, 0, -1, 0, 0, 0,},// gTbb|gNab >= 100*gTaa|gNab
            //{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, -1, 0,},// gTab|gNbb >= 100*gTaa|gNbb
            }
        };

        if (estAlpha) {
            return As[0];
        } else {
            return As[1];
        }

    }

    private int[][] getb(int[][] A) {
        int[][] b = new int[A.length][1];
        for (int i = 0; i < b.length; i++) {
            if (estAlpha && i == 1) {
                b[i] = new int[]{1};
            } else {
                b[i] = new int[]{0};
            }
        }

        return b;
    }

    private double[][] getQ() {
        double[][] q = {
            {0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1},};

        return q;
    }

    private double getLogLikelihood(double[] V) {
        double lr = 0;

        for (int i = 0; i < sc; i++) {
            double t = getLikelihoodAt(V, i);
            lr += Math.log(t);
        }

        return lr;
    }

    /*
     * private double[] getGradient(double[] V) { double[] g = new double[13];
     * Double stepSize = 1e-3;
     *
     *
     * if (estAlpha) { for (int i = 0; i < sc; i++) { double t = T(V, i, 1);
     * double s = 0;
     *
     * double[] denoTable = getDenoTable(af[i], V[0], minBnum, minBratio);
     * double[] denoTable2; if (V[0] + stepSize < 1) { denoTable2 =
     * getDenoTable(af[i], V[0] + stepSize, minBnum, minBratio); } else {
     * denoTable2 = getDenoTable(af[i], V[0] - stepSize, minBnum, minBratio); }
     *
     * for (int k = 0; k < 3; k++) { if (V[k + 1] == 0) { continue; }
     *
     * double u = 0; for (int j = 0; j < f[0][i].length; j++) { u +=
     * Math.log(f[0][i][j][k]); }
     *
     * for (int l = 0; l < 3; l++) { double x = u; double y = 0; // for (int j =
     * 0; j < f[0][i].length; j++) { // x += Math.log(f[0][i][j][k]); // } for
     * (int j = 0; j < f[1][i].length; j++) { double fk = f[1][i][j][k]; double
     * fl = f[1][i][j][l]; double h = h(V[0], fk, fl, k, l); x += Math.log(h); y
     * += (fk - fl) / h; }
     *
     * double deno = denoTable[k * 3 + l]; double deno2 = denoTable2[k * 3 + l];
     * double derideno = 0; if (V[0] + stepSize < 1) { derideno = (deno2 - deno)
     * / stepSize; } else { derideno = (deno - deno2) / stepSize; } //
     * if(derideno != 0) // System.out.println(y + "*" + derideno + "*" + deno +
     * "*" + (derideno + deno * y)); y = (derideno + deno * y); //s -=
     * getLogDenominator(k, l, V[0], f[1][i].length, i, minBnum); double z =
     * Math.pow(Math.E, x); // System.out.println(z + "\t" + deno); if (k == 0)
     * { //if (!(plotForMatlab && l == 1)) { g[4 + 3 * k + l] += (z * V[k + 1] /
     * t * deno);
     *
     * //} }
     *
     * if (k == l) { continue; } if (V[k * 3 + l + 4] > 0) { s += z * V[k + 1] *
     * V[k * 3 + l + 4] * y; } } } //if (!plotForMatlab) { //if(s<0)
     * System.out.println("YES"); g[0] += s / t; //}
     *
     * }
     * } else { for (int i = 0; i < sc; i++) { double t = T(V, i, 1); for (int k
     * = 0; k < 3; k++) { double gk = 0; for (int l = 0; l < 3; l++) { double x
     * = Math.pow(Math.E, logx[i][3 * k + l]); // float deno =
     * (float)Math.pow(Math.E, -getLogDenominator(k, l, V[0], f[1][i].length, i,
     * minBnum)); // x *= deno;
     *
     * if (V[k * 3 + l + 4] > 0) { gk += x * V[k * 3 + l + 4]; } if (V[k + 1] >
     * 0) { g[4 + 3 * k + l] += (x * V[k + 1] / t); }
     *
     * }
     * g[1 + k] += gk / t; } } }
     *
     *
     * double max = 0;
     *
     * for (int i = 0; i < g.length; i++) { max = Math.max(max, Math.abs(g[i]));
     * // if(Double.isNaN(g[i])) System.out.println("GShit"); }
     *
     * if (max > 0) { for (int i = 0; i < g.length; i++) { g[i] /= max; } }
     *
     * // MC.printVector(g); return g; }
     */
    private double[] getGradientAt(double[] V, int i, int version) {
        double[] g = new double[13];
        double stepSize = V[0]/100;
        if (estAlpha) {
            double t = getLikelihoodAt(V, i);
            double s = 0;

            double[] denoTable = getDenoTable(af[i], V[0], minBratio, maxBratio);
            double[] denoTable2;
            if (V[0] + stepSize < 1) {
                denoTable2 = getDenoTable(af[i], V[0] + stepSize, minBratio, maxBratio);
            } else {
                denoTable2 = getDenoTable(af[i], V[0] - stepSize, minBratio, maxBratio);
            }

            for (int k = 0; k < 3; k++) {
                if (V[k + 1] == 0) {
                    continue;
                }

                if (version == 2 && k != 0) {
                    break;
                }

                double u = 0;
                for (int j = 0; j < f[0][i].length; j++) {
                    u += Math.log(f[0][i][j][k]);
                }

                for (int l = 0; l < 3; l++) {

                    if (version == 0 && k == l) {
                        continue;
                    }


                    double x = u;
                    double y = 0;

                    for (int j = 0; j < f[1][i].length; j++) {
                        double fk = f[1][i][j][k];
                        double fl = f[1][i][j][l];
                        double h = h(V[0], fk, fl, k, l);
                        x += Math.log(h);
                        y += (fk - fl) / h;
                    }
                    double deno = denoTable[k * 3 + l];
                    if (version == 0) {
                        double deno2 = denoTable2[k * 3 + l];
                        double derideno = 0;
                        if (V[0] + stepSize < 1) {
                            derideno = (deno2 - deno) / stepSize;
                        } else {
                            derideno = (deno - deno2) / stepSize;
                        }

                        //if(version==0) 
                        y = (derideno + deno * y);
                    }

                    double z = Math.pow(Math.E, x);

                    if (version == 2 && k == 0) {
                        g[4 + 3 * k + l] += (z * V[k + 1] / t * deno);
                    }

                    if (version == 0 && V[k * 3 + l + 4] > 0) {
                        s += z * V[k + 1] * V[k * 3 + l + 4] * y;
                    }
                }
            }
            g[0] += s / t;
            //g[6] = 0;
        } else {
            double t = getLikelihoodAt(V, i);
            for (int k = 0; k < 3; k++) {
                if (version > 1 && version - 2 != k) {
                    continue;
                }

                double gk = 0;

                for (int l = 0; l < 3; l++) {
                    double x = Math.pow(Math.E, logx[i][3 * k + l]);

                    if (version == 1 && V[k * 3 + l + 4] > 0) {
                        gk += x * V[k * 3 + l + 4];
                    }
                    if (version != 1 && V[k + 1] > 0) {
                        g[4 + 3 * k + l] += (x * V[k + 1] / t);
                    }

                }
                if (version == 1) {
                    g[1 + k] += gk / t;
                }
            }
        }
        if (plotForMatlab) {
            g[0] = g[5] = 0;
        }
        return g;
    }

    private double[] makeZeroAccoringToVersion(double[] v, int version) {
        HashSet<Integer> nonZeroIndexes = new HashSet<Integer>();
        double[] d = new double[v.length];
        System.arraycopy(v, 0, d, 0, d.length);

        if (version == 0) {
            nonZeroIndexes.add(0);
        } else if (version == 1) {
            nonZeroIndexes.add(1);
            nonZeroIndexes.add(2);
            nonZeroIndexes.add(3);
        } else if (version == 2) {
            nonZeroIndexes.add(4);
            nonZeroIndexes.add(5);
            nonZeroIndexes.add(6);
        } else if (version == 3) {
            nonZeroIndexes.add(7);
            nonZeroIndexes.add(8);
            nonZeroIndexes.add(9);
        } else if (version == 4) {
            nonZeroIndexes.add(10);
            nonZeroIndexes.add(11);
            nonZeroIndexes.add(12);
        }

        for (int i = 0; i < d.length; i++) {
            if (!nonZeroIndexes.contains(i)) {
                d[i] = 0;
            }
        }

        if (estAlpha) {
            d[6] = 0;
        } else {
            d[0] = 0;
        }

        if (plotForMatlab) {
            d[0] = d[5] = 0;
        }

        return d;
    }

    private double[] getGradient(double[] V, int version) {
        double[] g = new double[13];

        if (estAlpha) {
            if (version != 0 && version != 2) {
                return null;
            }
        } else {
            if (version == 0) {
                return null;
            }
        }

        for (int i = 0; i < sc; i++) {
            double[] gi = getGradientAt(V, i, version);
            for (int n = 0; n < g.length; n++) {
                g[n] += gi[n];
            }

        }


        double max = 0;
        for (int i = 0;
                i < g.length;
                i++) {
            max = Math.max(max, Math.abs(g[i]));
        }
        if (max > 0) {
            for (int i = 0; i < g.length; i++) {
                g[i] /= max;
            }
        }
        // MC.printVector(g);

        //g = makeZeroAccoringToVersion(g, version);

        boolean hasNonZero = false;
        for (double gi : g) {
            if (gi != 0.0) {
                hasNonZero = true;
                break;
            }
        }

        if (!hasNonZero) {
            g = null;
        }

        return g;
    }

    private double[] getFeasibleDirection(double[] g, double[][] A1, double[][] Q, int version) {
        double[] d = new double[Q[0].length];

        LinearObjectiveFunction lf = new LinearObjectiveFunction(g, 0);

        Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();

        for (int i = 0; i < d.length; i++) {
            double[] a = new double[d.length];
            a[i] = 1;
            constraints.add(new LinearConstraint(a, Relationship.LEQ, 1));
            constraints.add(new LinearConstraint(a, Relationship.GEQ, -1));
        }
        // for(double z : g){
        //     System.out.println(z);
        // }

        // constraints.add(new LinearConstraint(g, Relationship.LEQ, 1));

        for (double[] a : Q) {
            constraints.add(new LinearConstraint(a, Relationship.EQ, 0));
        }

        for (double[] a : A1) {
            constraints.add(new LinearConstraint(a, Relationship.LEQ, 0));
        }
        PointValuePair solution = new SimplexSolver().optimize(lf, constraints, GoalType.MAXIMIZE, false);
        System.arraycopy(solution.getPoint(), 0, d, 0, d.length);
        d = makeZeroAccoringToVersion(d, version);
        //  MC.printVector(g);
        //  MC.printVector(d);

        return d;
    }

    private double getLambdaMax(double[] V, double[] d, double[][] A2, double[][] b2) {
        double m = 100;

        double[] bhat = MC.multiply(A2, V);
        double[] dhat = MC.multiply(A2, d);

        for (int i = 0; i < bhat.length; i++) {
            bhat[i] = b2[i][0] - bhat[i];
        }

        for (int i = 0; i < dhat.length; i++) {
            if (dhat[i] > 0) {
                m = Math.min(m, bhat[i] / dhat[i]);
            }
        }

        return m;
    }

    private double getLambda(double[] V, double[] d, double lambdaMax, double maxLR, int v) {
        double lambda;
        double optimalLambda = 0;
        double lr = -Double.MAX_VALUE;
        double[] tV;// = new double[V.length];

        int in = 0;

        double prevLR = maxLR;

        /*
         * while (lambda <= lambdaMax && in++ < 10) { tV = new float[V.length];
         *
         * for (int j = 0; j < tV.length; j++) { // if (d[j] != 0) { tV[j] =
         * Math.max(0, V[j] + (lambda * (float) d[j]));
         *
         * if(tV[j] <= delta){ tV[j] = 0; } }
         *
         * tV = normalizeV(tV); lr = getLogLikelihood(tV);
         *
         * if (prevLR < lr) { optimalLambda = lambda; prevLR = lr; }
         *
         * lambda *= 10; }
         */
        lambda = lambdaMax;


        while (in++ < 100) {
            tV = new double[V.length];

            for (int j = 0; j < tV.length; j++) {
                // if (d[j] != 0) {            
                tV[j] = (V[j] + (lambda * d[j]));
                if (tV[j] <= delta) {
                    tV[j] = 0;
                }

            }

            tV = normalizeV(tV);
            lr = getLogLikelihood(tV);

            if (prevLR <= lr) {
                optimalLambda = lambda;
                prevLR = lr;
                //break;
            }

            lambda *= .95f;
        }

//
//        //if (optimalLambda == 0) 
//        {
//            double p = optimalLambda*0.8, n = Math.max(optimalLambda*1.2, lambdaMax);
//            
//            
//            in = 0;
//            Random r = new Random();
//            while (p + delta < n && in++ < 20) {
//                optimalLambda = (p + n) / (2.0+r.nextDouble()*0.01);
//               
//                tV = new double[V.length];
//
//                for (int j = 0; j < tV.length; j++) {
//                    tV[j] = (V[j] + (optimalLambda * d[j]));
//
//                    if (tV[j] <= delta) {
//                        tV[j] = 0;
//                    }
//                }
//
//                tV = normalizeV(tV);
//
//                double[] g = getGradient(tV, v);
//                
//                double m = 0;
//                for (int j = 0; j< d.length; j++) {
//                    m += d[j] * g[j];
//                }
//
//                if (m > 0) {
//                    p = (4.0*p + optimalLambda) / 5.0;
//                } else if (m == 0) {
//                    break;
//                } else if (m < 0) {
//                    n = (optimalLambda + 4.0*n) / 5.0;                   
//                }
//            }
//        }



        return optimalLambda;
    }

    private double feasibleDirection(double[] V, int[][] A, int[][] b, double[][] Q) {
        double[] lambda = new double[5];
        // float prevLR = -Float.MAX_VALUE;
        double lr = getLogLikelihood(V);
        double lambdaMax = 1;
        // System.out.println("before " + lr);
        int in = 0;

        while (in++ < 5) {// lambda > Math.min(lambdaMax / 3, delta)  lambda > delta &&

            for (int v = 4; v >= 0; v--) {
                ArrayList<Integer> br = getBindingRaws(V, A, b);
                double[] g = getGradient(V, v);
                if (g == null) {
                    continue;
                }
                double[] d = getFeasibleDirection(g, getA1(A, br), Q, v);

                lambdaMax = getLambdaMax(V, d, getAb2(A, br), getAb2(b, br));
                lambda[v] = getLambda(V, d, lambdaMax, lr, v);

                for (int j = 0; j < d.length; j++) {
                    V[j] = (V[j] + (lambda[v] * d[j]));//Math.min(Math.max(V[j] + (lambda * (float) d[j]), 1e-20f), 1 - 1e-20f);
                    if (V[j] <= delta) {
                        V[j] = 0;
                    }
                }

                V = normalizeV(V);
                lr = getLogLikelihood(V);

                /*
                 * float[] tV = new float[V.length]; for (int i = 0; i <
                 * V.length; i++) { tV[i] = V[i]; } // tV[0] = 0.9f; tV[4] =
                 * 0.85f; tV[5] = 0.15f; tV[6] = 0.00000001f; tV =
                 * normalizeV(tV); float tlr = getLogLikelihood(tV);
                 * System.out.println("\n" + v + "\t" + V[0] + "\t" + (lr > tlr)
                 * + "\t" + d[0] + "\t" + lambda[v] + "\t" + lambdaMax);
                 * MC.printVector(g); MC.printVector(d);
                 */
                // MC.printVector(V);
                //    System.out.println(v + "\t" + lambdaMax + "\t" + lambda[v] + "\t : " + lr);
                //    MC.printVector(V);
                //   System.out.println(v + "\t" + lambdaMax + "\t" + lambda[v]);
                //    MC.printVector(g);
                //    MC.printVector(d);

            }

            boolean converged = true;

            for (double lm : lambda) {
                if (lm != 0.0) {
                    converged = false;
                }
            }

            if (!plotForMatlab && converged) {
                break;
            }
        }

        return lr;
    }

    private double[] normalizeV(double[] V) {
        double ss = 0;
        for (int i = 0; i < 3; i++) {
            V[i + 1] = Math.max(V[i + 1], 0);
            ss += V[i + 1];
        }

        for (int i = 0; i < 3; i++) {
            V[i + 1] /= ss;
        }

        for (int i = 0; i < 3; i++) {
            double s = 0;
            for (int j = 0; j < 3; j++) {
                V[i * 3 + j + 4] = Math.max(V[i * 3 + j + 4], 0);
                s += V[i * 3 + j + 4];
            }
            for (int j = 0; j < 3; j++) {
                if (!plotForMatlab) {
                    V[i * 3 + j + 4] /= s;
                }

            }
        }

        if (plotForMatlab) {
            V[6] = V[4] / 1000000f;
            V[4] = 1 - V[5] - V[6];
        }

        V[0] = Math.min(V[0], 1 - 1e-2);
        return V;
    }

    private double[] getInitialV(double seed) {
        double[] newV;

        double[] G1 = {
            5e-1, 5e-1, 1e-7,
            1e-40, 1e-40, 1e-40,
            1e-40, 1e-40, 1e-40,};

        double[] G2 = {
            1, 5e-5, 1e-9,
            1e-20, 5e-4, 1e-9,
            1e-20, 1e-20, 1e-10};


        if (estAlpha) {
            // if (this.V == null) {
            newV = getV(seed, G1);
            //} else {
            //    newV = new float[this.V.length];
            //    System.arraycopy(this.V, 0, newV, 0, newV.length);
            //V[0] = alpha;
            //}
        } else {
            // if (this.V == null) {
            newV = getV(getInitialAlpha(), G2);
            // } else {
            //     newV = new float[this.V.length];
            //      System.arraycopy(this.V, 0, newV, 0, newV.length);
            // }
        }
        // System.out.println("**");
        // MC.printVector(newV);
        //3.97e-01 1.00e+00 1.11e-20 1.11e-20 4.44e-01 5.56e-01 1.11e-05 3.33e-01 3.33e-01 3.33e-01 3.33e-01 3.33e-01 3.33e-01 
        //

        double[] tV = new double[newV.length];
        System.arraycopy(newV, 0, tV, 0, newV.length);

        int cntr = 100;
        while (true) {
            Random r = new Random();
            for (int i = 1; i < newV.length; i++) {
                double t = newV[i];
                double add = (r.nextDouble() - 0.5f) * Math.min(0.1f, Math.abs(t));
                t += add;

                //if(i==0)System.out.println(tV);

                newV[i] = Math.min(Math.max(t, 1e-20f), 1 - 1e-20f);
            }
            // if (estAlpha) {
            //     newV[0] += (r.nextFloat() - 0.5f) * getInitialAlpha() * 0.5f;
            //     newV[0] = Math.min(Math.max(newV[0], 1e-20f), 1 - 1e-20f);
            // }

            if (plotForMatlab) {
                newV[0] = this.V[0];
                newV[5] = this.V[5];

            }
            newV = normalizeV(newV);

            if (getBindingRaws(newV, getA(), this.getb(getA())).isEmpty() || cntr-- < 0) {
                break;
            } else {
                System.arraycopy(tV, 0, newV, 0, tV.length);
            }
        }
        return newV;
    }

    private double[] getV(double alpha, double[] G) {
        double[] tV = new double[13];

        tV[0] = alpha;

        for (int i = 0; i < 3; i++) {
            tV[i + 1] = 0;
            for (int j = 0; j < 3; j++) {
                tV[i + 1] += G[i * 3 + j];
            }
            tV[i + 1] = Math.min(Math.max(1e-20f, tV[i + 1]), 1 - 1e-2f);
        }


        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (tV[i + 1] > 0) {
                    tV[i * 3 + j + 4] = G[i * 3 + j] / tV[i + 1];
                }
                tV[i * 3 + j + 4] = Math.min(Math.max(1e-20f, tV[i * 3 + j + 4]), 1 - 1e-20f);
            }

        }

        return normalizeV(tV);
    }

    private void clear() {
        this.f = null;
        this.logx = null;
        this.af = null;
        this.sc = 0;
    }

    public double[] getEstimationResults() {
        return V;
    }

    public void setEstimationParameters(double[] V) {
        this.V = V;
    }

    public void feed(boolean[][] T, float[][] erT, float[][] emT, boolean[][] N, float[][] erN, float[][] emN) {


        //feedSingle(boolean[] T, float[] erT, float[] emT, boolean[] N, float[] erN, float[] emN) 
        for (int i = 0; i < T.length; i++) {
            feedSingle(T[i], erT[i], emT[i], N[i], erN[i], emN[i]);
        }

        //int minBnum = (int) (Math.ceil(af.length * minBratio)+0.01);

        System.out.println("Sample Number : " + sc);
        /*
         * if (estAlpha) { this.f = new double[2][T.length][][]; this.af = new
         * double[T.length][][];
         *
         * for (int i = 0; i < T.length; i++) { this.f[0][sc + i] = new
         * double[N[i].length][3]; for (int j = 0; j < N[i].length; j++) { for
         * (int k = 0; k < 3; k++) { this.f[0][sc + i][j][k] = f(k, erN[i][j],
         * emN[i][j], N[i][j], beta); } }
         *
         * this.f[1][sc + i] = new double[T[i].length][3]; this.af[sc + i] = new
         * double[T[i].length][3];
         *
         * double tt = 0; for (int j = 0; j < T[i].length; j++) { for (int k =
         * 0; k < 3; k++) { this.f[1][sc + i][j][k] = f(k, erT[i][j], emT[i][j],
         * T[i][j], beta);
         *
         *
         * if (T[i][j]) { this.af[sc + i][j][k] = this.f[1][sc + i][j][k]; }
         * else { this.af[sc + i][j][k] = 1 - this.f[1][sc + i][j][k]; tt++; } }
         * }
         *
         * // if(tt/3f/T[i].length<minBratio){ // System.out.println("**** " +
         * T[i].length +"\t" + minBratio); // int ss = 0; // for(boolean b :
         * T[i]){ // System.out.print(b+" "); // if(!b) ss++; // } //
         * System.out.println(ss); // } } } else { for (int i = 0; i < T.length;
         * i++) {
         *
         * double[][] fs = new double[3][T[i].length]; double[][] afs = new
         * double[T[i].length][3];
         *
         * for (int k = 0; k < 3; k++) { for (int j = 0; j < T[i].length; j++) {
         * fs[k][j] = f(k, erT[i][j], emT[i][j], T[i][j], beta); if (T[i][j]) {
         * afs[j][k] = fs[k][j]; } else { afs[j][k] = 1 - fs[k][j]; } } }
         *
         * //double[] denoTable = getDenoTable(afs, this.getInitialAlpha(),
         * this.minBnum, minBratio);
         *
         * for (int k = 0; k < 3; k++) { for (int l = 0; l < 3; l++) { int index
         * = 3 * k + l; for (int j = 0; j < T[i].length; j++) { double fk =
         * fs[k][j];//f(k, erT[i][j], emT[i][j], T[i][j]); double fl =
         * fs[l][j];//f(l, erT[i][j], emT[i][j], T[i][j]); logx[sc + i][index]
         * += Math.log(h(this.getInitialAlpha(), fk, fl, k, l)); }
         *
         * //logx[sc + i][index] += Math.log(denoTable[index]); } } } for (int
         * i = 0; i < N.length; i++) { double[][] fs = new
         * double[3][N[i].length]; for (int j = 0; j < N[i].length; j++) { for
         * (int k = 0; k < 3; k++) { fs[k][j] = f(k, erN[i][j], emN[i][j],
         * N[i][j], beta); } } for (int j = 0; j < N[i].length; j++) { for (int
         * k = 0; k < 3; k++) { double z = Math.log(fs[k][j]); for (int l = 0; l
         * < 3; l++) { logx[sc + i][3 * k + l] += z; } } } } } sc += T.length;
         */

    }

    /*
     * private ArrayList<Integer> getQualifiedIndices(float[] er, float
     * threshold){ ArrayList<Integer> indices = new ArrayList<Integer>();
     *
     * for(int i=0; i<er.length; i++){ if(er[i]<threshold) indices.add(i); }
     *
     * return indices; }
     */
    public void feedSingle(boolean[] T, float[] erT, float[] emT, boolean[] N, float[] erN, float[] emN) {

        if (estAlpha) {
            int minBnum = (int) (Math.ceil(T.length * minBratio) + 0.01);
            double cntr = 0;
            for (boolean tt : T) {
                if (!tt) {
                    cntr++;
                }
            }


            if (cntr < minBnum) {
               return;
            }
            
            if (maxBratio > 0) {              
                if (cntr / T.length > maxBratio) {
                    return;
                }
            }

           
            this.f[0][sc] = new double[N.length][3];
            for (int j = 0; j < N.length; j++) {
                for (int k = 0; k < 3; k++) {
                    this.f[0][sc][j][k] = f(k, erN[j], emN[j], N[j],
                            beta);
                }

            }
            this.f[1][sc] = new double[T.length][3];
            this.af[sc] = new double[T.length][3];
            for (int j = 0; j < T.length; j++) {
                for (int k = 0; k < 3; k++) {
                    this.f[1][sc][j][k] = f(k, erT[j], emT[j], T[j], beta);
                    if (T[j]) {
                        this.af[sc][j][k] = this.f[1][sc][j][k];
                    } else {
                        this.af[sc][j][k] = 1 - this.f[1][sc][j][k];
                    }
                }
            }
        } else {
            double[][] fs = new double[3][T.length];
            double[][] afs = new double[T.length][3];

            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < T.length; j++) {
                    fs[k][j] = f(k, erT[j], emT[j], T[j], beta);
                    if (T[j]) {
                        afs[j][k] = fs[k][j];
                    } else {
                        afs[j][k] = 1 - fs[k][j];
                    }
                }
            }

            //double[] denoTable = getDenoTable(afs, this.getInitialAlpha() , this.minBnum  , minBratio            );

            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    int index = 3 * k + l;
                    for (int j = 0; j < T.length; j++) {
                        double fk = fs[k][j];
                        double fl = fs[l][j];//f(l, erT[i][j], emT[i][j], T[i][j]
                        logx[sc][index] += Math.log(h(this.getInitialAlpha(), fk, fl, k, l));
                    }
                    //logx[sc][index] += Math.log(denoTable[index]);
                }
            }


            double[][] fsn = new double[3][N.length];
            for (int j = 0; j < N.length; j++) {
                for (int k = 0; k < 3; k++) {
                    fsn[k][j] = f(k, erN[j], emN[j], N[j], beta);
                }
            }
            for (int j = 0; j < N.length; j++) {
                for (int k = 0; k < 3; k++) {
                    double z = Math.log(fsn[k][j]);
                    for (int l = 0; l < 3; l++) {
                        logx[sc][3 * k + l] += z;
                    }
                }
            }

        }



        sc += 1;
    }

    public double estimate() {

        int[][] A = this.getA();
        int[][] b = this.getb(A);
        double[][] Q = this.getQ();
//        float maxLR = -Float.MAX_VALUE;
        double maxLR = Double.NEGATIVE_INFINITY;


//        System.out.print("Iteration .. ");
        int cntr = 10;
        int numIteration = (plotForMatlab ? 1 : (estAlpha ? 5 : 3));
        for (int k = 0; k < numIteration; k++) {//
            double[] tV = getInitialV((double) k / numIteration + .01);//

            double LR = this.feasibleDirection(tV, A, b, Q);
            //float mLR = LR;
            //float maxAlpha = tV[0];

            /*
             * if(estAlpha){ for(float a=0.01f;a<1;a+=0.01f){ tV[0] = a; LR =
             * this.getLogLikelihood(tV); if(mLR < LR){ mLR = LR; maxAlpha = a;
             * } }
             *
             * tV[0] = maxAlpha; }
             */

            //if (maxLR <= LR) {
                //System.out.print("%*");
            //} else {
                //System.out.print("%");
            //}
            //if (tV != null) {
                //MC.printVector(tV);
                //System.out.println("%Likelihood : " + LR);
            //}


            if (maxLR <= LR) {
                this.V = tV;
                maxLR = LR;
            }


            if (cntr-- <= 0) {
                break;
            }
        }


        /*
         * if (estAlpha && !plotForMatlab) { BufferedWriter alphaSpaceWriter =
         * null; try { alphaSpaceWriter = new BufferedWriter(new
         * FileWriter(Param.alphaSpaceFile)); float at = V[0]; float called =
         * getLogLikelihood(V); boolean error = false; for (float a = 0.01f; a <
         * 0.99f; a += 0.01f) { float[] tV = new float[V.length];
         * System.arraycopy(V, 0, tV, 0, tV.length); tV[0] = a; float other =
         * getLogLikelihood(tV); if (other > called) { error = true; }
         * alphaSpaceWriter.write(String.valueOf(a) + "\t" +
         * String.valueOf(other)); alphaSpaceWriter.write("\n"); if (V[0] > a &&
         * V[0] < a + 0.01f) { alphaSpaceWriter.write(String.valueOf(V[0]) +
         * "\t" + String.valueOf(called)); alphaSpaceWriter.write("\n"); } } if
         * (error) { System.out.println("Estimation error occured. (" +
         * String.valueOf(at) + ":" + String.valueOf(called) + ") see " +
         * Param.alphaSpaceFile.getName()); alphaSpaceWriter.close(); //
         * System.exit(-1); } } catch (IOException ex) {
         * Logger.getLogger(MLE.class.getName()).log(Level.SEVERE, null, ex); }
         * finally { try { alphaSpaceWriter.close(); } catch (IOException ex) {
         * Logger.getLogger(MLE.class.getName()).log(Level.SEVERE, null, ex); }
         * } }
         */
        /*
         * if (estAlpha) {
         *
         * // System.out.println("LR with estimated alpha " +
         * getLogLikelihood(V)); // Param.report.write("LR with estimated alpha
         * " + String.valueOf(getLogLikelihood(V))); for (float a = V[0] / 10; a
         * < Math.min(V[0] * 2, 0.99f); a += V[0] / 10) { float[] tV = new
         * float[V.length];
         *
         * System.arraycopy(V, 0, tV, 0, tV.length); //tV[0] = 0.01f; //
         * System.out.println("LR with real alpha " + getLogLikelihood(tV));
         * tV[0] = a; // System.out.println("alpha : " + a + "\tLR : " +
         * getLogLikelihood(tV)); } }
         */

        for (int i = 0; i < V.length; i++) {
            V[i] = Math.max(V[i], Double.MIN_VALUE * 2);
            V[i] = Math.min(V[i], 1);
        }
        double sumk = 0;
        double[] suml = new double[3];
        for (int k = 0;
                k < 3; k++) {
            sumk += V[k + 1];
            for (int l = 0; l < 3; l++) {
                suml[k] += V[4 + k * 3 + l];
            }
        }
        for (int k = 0;
                k < 3; k++) {
            V[k + 1] /= sumk;
            for (int l = 0; l < 3; l++) {
                V[4 + k * 3 + l] /= suml[k];
            }
        }

        if (estAlpha) {
            System.out.println("%alpha estimation variance : " + getCovarianceMatrix(V)[0][0]);
        }

        if (!plotForMatlab) {
            clear();
        }


        return maxLR;
    }

    /*
     * private float[] getNormalPosteriorGenotypeDist(boolean[] T, float[] erT,
     * float[] emT, boolean[] N, float[] erN, float[] emN, float[] V) { double[]
     * tV = new double[V.length];
     *
     * for(int i=0;i<tV.length;i++){ tV[i] = V[i]; } // System.arraycopy(V, 0,
     * tV, 0, tV.length);
     *
     * double maxExp0 = Double.MIN_EXPONENT;
     *
     * for (int k = 0; k < 3; k++) { double prod = 0; for (int j = 0; j <
     * N.length; j++) { prod += Math.log(f(k, erN[j], emN[j], N[j])); }
     *
     * //double t = 0;
     *
     *
     * for (int l = 0; l < 3; l++) { double prod2 = 0; for (int j = 0; j <
     * T.length; j++) { float fk = f(k, erT[j], emT[j], T[j]); float fl = f(l,
     * erT[j], emT[j], T[j]); prod2 += Math.log(h(V[0], fk, fl, k, l)); } double
     * c = prod2 + Math.log(V[4 + k * 3 + l]);
     *
     * tV[4 + k * 3 + l] = c;
     *
     * }
     *
     * double t = 0;
     *
     * for (int l = 0; l < 3; l++) { tV[4 + k * 3 + l] = Math.pow(Math.E, tV[4 +
     * k * 3 + l] + 100); t += tV[4 + k * 3 + l]; }
     *
     * for (int l = 0; l < 3; l++) { tV[4 + k * 3 + l] /= t; }
     *
     *
     * double c = Math.log(t) + Math.log(V[1 + k]) + prod;
     *
     * maxExp0 = Math.max(maxExp0, c);
     *
     * tV[k + 1] = c;
     *
     * }
     *
     * double t = 0;
     *
     * for (int k = 0; k < 3; k++) { tV[k + 1] = Math.pow(Math.E, tV[k + 1] -
     * maxExp0); t += tV[k + 1]; } for (int k = 0; k < 3; k++) { tV[k + 1] /= t;
     * }
     *
     * float[] rV = new float[tV.length]; for(int i = 0; i<tV.length; i++){
     * rV[i] = (float)tV[i]; } * return rV; }
     *
     * private float[] getTumorPosteriorGenotypeDist(boolean[] T, float[] erT,
     * float[] emT, float[] V) {
     *
     * float[] tV = new float[V.length];
     *
     * System.arraycopy(V, 0, tV, 0, tV.length);
     *
     * for (int k = 0; k < 3; k++) { float sum = 0;
     *
     * for (int l = 0; l < 3; l++) { double prod = 0;
     *
     * for (int j = 0; j < T.length; j++) { float fk = f(k, erT[j], emT[j],
     * T[j]); float fl = f(l, erT[j], emT[j], T[j]); prod += Math.log(h(V[0],
     * fk, fl, k, l)); }
     *
     * float c = Math.max(Float.MIN_VALUE, (float)(Math.pow(Math.E, prod) * V[4
     * + k * 3 + l]));
     *
     * tV[4 + k * 3 + l] = c;// * V[1 + k]; sum += c; } for (int l = 0; l < 3;
     * l++) { tV[4 + k * 3 + l] /= sum; } }
     *
     * return tV; }
     */
    private float[][] getGenoTypeMatrix(double[] V) {
        float[][] G = new float[3][3];
        float s = 0;

        V = normalizeV(V);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                G[i][j] = (float) (V[i * 3 + j + 4] * V[i + 1]);
                G[i][j] = Math.min(1, Math.max(0, G[i][j]));
                s += G[i][j];
            }
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                G[i][j] /= s;
            }
        }
        return G;
    }

    public float[][] predictGenotype2(boolean[] T, float[] erT, float[] emT, boolean[] N, float[] erN, float[] emN) {

        double[] tV = new double[V.length];
        //ouble beta = 0.99f;
        for (int i = 0; i < tV.length; i++) {
            tV[i] = V[i];
        }

        double maxLR = -Double.MAX_VALUE;
        for (int b = 0; b <= 5; b++) {
            double lr = 0;
            for (int k = 0; k < 3; k++) {
                double prod = 20;
                for (int j = 0; j < N.length; j++) {
                    prod += Math.log(f(k, erN[j], emN[j], N[j], .2f * b));
                }
                lr += Math.pow(Math.E, prod) * tV[k + 1];
            }
            if (maxLR < lr) {
                maxLR = lr;
                beta = .2 * b;
            }
        }

        beta = Math.max(1e-10, Math.min(beta, 1 - 1e-7));

        /*
         * if(beta < 0.3){ System.out.print(beta + "\t"); for(boolean n : N)
         * System.out.print(n+" "); System.out.println(); }
         */

        double maxExp0 = Double.MIN_EXPONENT;

        double[][] fs = new double[3][T.length];
        //float[][] afs = new float[T.length][3];

        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < T.length; j++) {
                fs[k][j] = f(k, erT[j], emT[j], T[j], beta);
                /*
                 * if (T[j]) { afs[j][k] = fs[k][j]; } else { afs[j][k] = 1 -
                 * fs[k][j]; }
                 */
            }
        }

        //float[] denoTable = getDenoTable(afs, V[0], this.minBnum, minBratio);
        //MC.printVector(denoTable);
        for (int k = 0; k < 3; k++) {
            double prod = 0;
            for (int j = 0; j < N.length; j++) {
                prod += Math.log(f(k, erN[j], emN[j], N[j], beta));
            }

            double t = 0;

            for (int l = 0; l < 3; l++) {
                double prod2 = 100;
                for (int j = 0; j < T.length; j++) {
                    double fk = fs[k][j];//f(k, erT[j], emT[j], T[j]);
                    double fl = fs[l][j];// f(l, erT[j], emT[j], T[j]);
                    prod2 += Math.log(h(V[0], fk, fl, k, l));
                }

                tV[4 + k * 3 + l] = Math.pow(Math.E, prod2) * V[4 + k * 3 + l];// * denoTable[k * 3 + l];
                t += tV[4 + k * 3 + l];
            }

            for (int l = 0; l < 3; l++) {
                tV[4 + k * 3 + l] /= t;
            }

            double c = Math.log(t) + Math.log(V[1 + k]) + prod;

            maxExp0 = Math.max(maxExp0, c);

            tV[k + 1] = c;

        }

        double t = 0;

        for (int k = 0; k < 3; k++) {
            tV[k + 1] = Math.pow(Math.E, tV[k + 1] - maxExp0);
            t += tV[k + 1];
        }
        for (int k = 0; k < 3; k++) {
            tV[k + 1] /= t;
        }

        double[] rV = new double[tV.length];
        for (int i = 0; i < tV.length; i++) {
            rV[i] = tV[i];
        }

        float[][] g = getGenoTypeMatrix(rV);

        return g;
    }

    public float[][] predictGenotype3(boolean[] T, float[] erT, float[] emT, boolean[] N, float[] erN, float[] emN) {

        double[] tV = new double[V.length];
        //System.arraycopy(V, 0, tV, 0, tV.length);

        double maxLR = -Double.MAX_VALUE;
        // double factor = 1;

        //for (int fa = 1; fa <= 50; fa += 2) {
        for (int b = 0; b <= 5; b++) {
            double lr = 0;
            for (int k = 0; k < 3; k++) {
                double prod = 20;
                for (int j = 0; j < N.length; j++) {
                    prod += Math.log(f(k, erN[j], emN[j], N[j], .2f * b));
                }
                lr += Math.pow(Math.E, prod) * V[k + 1];
            }
            if (maxLR < lr) {
                maxLR = lr;
                beta = .2 * b;
                // factor = fa;
            }
        }
        // }
        beta = Math.max(1e-10, Math.min(beta, 1 - 1e-7));


        double[][] ft = new double[3][T.length];

        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < T.length; j++) {
                ft[k][j] = f(k, erT[j], emT[j], T[j], beta);
            }
        }

        double[][] fn = new double[3][N.length];

        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < N.length; j++) {
                fn[k][j] = f(k, erN[j], emN[j], N[j], beta);
            }
        }

        double[][] ll = new double[3][3];
        double[][] hgg = new double[3][3];

        for (int k = 0; k < 3; k++) {
            double u = 100.0;
            for (int j = 0; j < N.length; j++) {
                u += Math.log(fn[k][j]);
            }

            tV[k + 1] = Math.pow(Math.E, u) * V[k + 1];

            for (int l = 0; l < 3; l++) {
                double s = u;

                for (int j = 0; j < T.length; j++) {
                    double fk = ft[k][j];
                    double fl = ft[l][j];
                    s += Math.log(h(V[0], fk, fl, k, l));
                }

                ll[k][l] = Math.pow(Math.E, s);
                hgg[k][l] = ll[k][l] * V[4 + k * 3 + l];

            }
        }

        /*
         * double t = 0; double[] s = new double[3];
         *
         * for (int k = 0; k < 3; k++) { for (int l = 0; l < 3; l++) { t +=
         * tgg[k][l]; s[k] += hgg[k][l]; } }
         */

        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                //tgg[k][l]/=t;
                //tV[k + 1] += tgg[k][l];
                // hgg[k][l]/=s[k];
                tV[4 + k * 3 + l] = hgg[k][l];
            }
        }

        float[][] g = getGenoTypeMatrix(tV);

        return g;
    }

    public float[][] predictGenotype(boolean[] T, float[] erT, float[] emT, boolean[] N, float[] erN, float[] emN) {

        /*
         * float nb = 0, lowq = 0; for (int j = 0; j < N.length; j++) { if
         * (!N[j]) { nb++; if (emN[j] > 0.01) { lowq++; } } }
         *
         * if (lowq / nb > .2) { return null; }
         */
        double[] tV = new double[V.length];

        //double[] factors = new double[3];
        //double factor = 1;

        //for (double fa = 1; fa < 1e3; fa *= 2) {
        double maxLR = -Double.MAX_VALUE;
        for (int b = 0; b <= 5; b++) {
            double lr = 0;
            for (int k = 0; k < 3; k++) {
                double prod = 20;
                for (int j = 0; j < N.length; j++) {
                    prod += Math.log(f(k, erN[j], emN[j], N[j], .2f * b));
                }
                lr += Math.pow(Math.E, prod) * V[k + 1];
            }
            if (maxLR < lr) {
                maxLR = lr;
                beta = .2 * b;
            }
        }

        beta = Math.max(1e-10, Math.min(beta, 1 - 1e-7));
        // }
        //beta = Math.max(1e-10, Math.min(beta, 1 - 1e-7));


        double[][] ft = new double[T.length][3];
        double[][] aft = new double[T.length][3];

        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < T.length; j++) {
                ft[j][k] = f(k, erT[j], emT[j], T[j], beta);
                if (T[j]) {
                    aft[j][k] = ft[j][k];
                } else {
                    aft[j][k] = 1 - ft[j][k];
                }
            }
        }

        double[][] fn = new double[N.length][3];
        double[][] afn = new double[N.length][3];

        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < N.length; j++) {
                fn[j][k] = f(k, erN[j], emN[j], N[j], beta);
                if (N[j]) {
                    afn[j][k] = fn[j][k];
                } else {
                    afn[j][k] = 1 - fn[j][k];
                }
            }
        }

        //double[][] ll = new double[3][3];
        double[][] tgg = new double[3][3];
        double[][] hgg = new double[3][3];
        //double[][] maxll = new double[3][3];
        double[] minU = new double[3];
        double[][] minS = new double[3][3];

        for (int k = 0; k < 3; k++) {
            double u = 100.0;
            for (int j = 0; j < N.length; j++) {
                double tu = Math.log(fn[j][k]);
                u += tu;
                minU[k] += -tu;
            }

            for (int l = 0; l < 3; l++) {
                double s = u;

                for (int j = 0; j < T.length; j++) {
                    double fk = ft[j][k];
                    double fl = ft[j][l];
                    double ts = Math.log(h(V[0], fk, fl, k, l));
                    s += ts;
                    minS[k][l] += -ts;
                }


                double lr = Math.pow(Math.E, s);
                hgg[k][l] = lr * V[4 + k * 3 + l];
                tgg[k][l] = lr * V[4 + k * 3 + l] * V[k + 1];

                hgg[k][l] = Math.max(hgg[k][l], Double.MIN_VALUE);
                tgg[k][l] = Math.max(tgg[k][l], Double.MIN_VALUE);
            }
        }



        /*
         * double t = 0; double[] s = new double[3];
         *
         * for (int k = 0; k < 3; k++) { for (int l = 0; l < 3; l++) { t +=
         * tgg[k][l]; s[k] += hgg[k][l]; } }
         */

        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                //tgg[k][l]/=t;
                tV[k + 1] += tgg[k][l];
                // hgg[k][l]/=s[k];
                tV[4 + k * 3 + l] = hgg[k][l];
            }
        }

        float[][] g = getGenoTypeMatrix(tV);

        double maxG = 0;
        int kk = 0, ll = 0;
        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                if (maxG < g[k][l]) {
                    maxG = g[k][l];
                    kk = k;
                    ll = l;
                }
            }
        }

        // System.out.println(minU[kk]+"\t" + minS[kk][ll]);
        double pval = Math.max(1 - getPValue(afn, 0, minU[kk], kk, kk) * getPValue(aft, V[0], minS[kk][ll], kk, ll), 1e-100);
//     

        // if(tt<Math.log(.3)){
        //      MC.printMatrix(g);
        //      System.out.println(maxll/(T.length + N.length) + "\t" + Math.log(.5));
        //  }
        beta = Math.log10(pval); // the higher the better..
        //if(tt<Math.log(.3)){// g[0][0] + g[0][1] + g[0][2] > 0.5 && (nb >= 3 && nb/(float)N.length >= 0.05)){
        // System.out.println("Filtered " + tt + "\t" + Math.log(.3));
        //   return null;
        // }
        //  if (beta < -2) {
        //       return null;
        //  }

        return g;
    }

    public double getBeta() {
        return beta;
    }

    public double getEstimatedAlpha() {
        return V[0];
    }

    public float[][] getEstimatedTransitionMatrix() {
        return getGenoTypeMatrix(V);
    }

    public void printResult() {
        System.out.println("Estimated alpha : " + getEstimatedAlpha());
        System.out.println("Estimated transition matrix :");
        MC.printMatrix(getEstimatedTransitionMatrix());


        //   V[0] = alpha;
        //  System.out.println("LR with true alpha  : " + getLogLikeliHood(readData, errorRates, V));
        //for(float v : V)
        //    System.out.println(v);
        //  System.out.println("LR with estimated alpha " + getLikeliHood(readData, errorRates, V));
        //   V[0] = alpha;
        //  System.out.println("LR with true alpha  : " + getLogLikeliHood(readData, errorRates, V));


    }

    public void printLikelihoodsForMatlab() {
        plotForMatlab = true;

        for (float alpha = 0.0f; alpha <= 1.01f; alpha += 0.05f) {
            for (float gab = 0.0f; gab <= 1.01f; gab += 0.05f) {
                this.V = new double[13];
                this.V[0] = Math.max(Math.min(0.99f, alpha), 0.001f);
                this.V[5] = Math.max(Math.min(0.99f, gab), 0.001f);
                System.out.println(alpha + "\t" + gab + "\t" + this.estimate());
            }
        }

        clear();

        plotForMatlab = false;
    }

    public double x(int read_length, int edit_distance, float error, float germ, float som) {
        float p = germ + som + error;
        double x1 = 0.0f;
        for (int i = 0; i <= edit_distance; i++) {
            int coef = (int) ArithmeticUtils.binomialCoefficient(read_length - 1, i);
            x1 += coef * Math.pow(p, i) * Math.pow(1 - p, read_length - i - 1);
        }
        return 1 - x1;
    }

    public double y(int read_length, int edit_distance, float error, float germ, float som) {
        float p = germ + som + error;
        double y1 = 0.0f;
        for (int i = 0; i <= edit_distance - 1; i++) {
            int coef = (int) ArithmeticUtils.binomialCoefficient(read_length - 1, i);
            y1 += coef * Math.pow(p, i) * Math.pow(1 - p, read_length - i - 1);
        }
        return 1 - y1;
    }

    public void setxy(int read_length, int edit_distance, float error, float germ, float som) {
        this.x = x(read_length, edit_distance, error, germ, som);
        this.y = y(read_length, edit_distance, error, germ, som);
    }

    /**
     * @return the alpha
     */
    public double getInitialAlpha() {
        return alpha;
    }

    /**
     * @param alpha the alpha to set
     */
    public void setAlpha(double alpha) {
        V[0] = alpha;
    }
}
