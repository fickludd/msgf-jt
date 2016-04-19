/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package virmid;

/**
 *
 * @author rachmani
 */
public class Pile {

    private String chrom;
    private int pos;
    private String ref;
    private int depth;
    private String call;
    private String qual;
    private String key;
    private String mappingqual;
    private boolean[] maskv;
    private char[] charactorv;
    public char major;
    public char minor;
    private int modifiedLength;
    private boolean valid = false;
    private boolean hasb = false;
    private String raw;

    public Pile() {
        init();
    }

    public Pile(String raw) {
        init();
        String[] fields = raw.split("\t");
        if (fields.length >= 6) {
            this.chrom = fields[0];
            this.pos = Integer.valueOf(fields[1]);
            this.ref = fields[2];
            this.depth = Integer.valueOf(fields[3]);
            this.call = fields[4];
            this.qual = fields[5];
            this.key = this.chrom + ":" + String.valueOf(this.pos);
            this.raw = raw;
            if (fields.length >= 7) {
                this.mappingqual = fields[6];
            } else {
                this.mappingqual = null;
            }
            this.getCharactorVector();
            this.eliminateMimorAlleles();
            this.valid = true;
        }
    }

    private void init() {
        this.setChrom("");
        this.setPos(0);
        this.setRef("");
        this.setDepth(0);
        this.setCall("");
        this.setQual("");

    }

    /**
     * @return the chrom
     */
    public String getChrom() {
        return chrom;
    }

    /**
     * @param chrom the chrom to set
     */
    public void setChrom(String chrom) {
        this.chrom = chrom;
    }

    /**
     * @return the pos
     */
    public int getPos() {
        return pos;
    }

    /**
     * @param pos the pos to set
     */
    public void setPos(int pos) {
        this.pos = pos;
    }

    /**
     * @return the ref
     */
    public String getRef() {
        return ref;
    }

    /**
     * @param ref the ref to set
     */
    public void setRef(String ref) {
        this.ref = ref;
    }

    /**
     * @return the depth
     */
    public int getDepth() {
        return depth;
    }

    /**
     * @param depth the depth to set
     */
    public void setDepth(int depth) {
        this.depth = depth;
    }

    /**
     * @return the call
     */
    public String getCall() {
        return call;
    }

    /**
     * @param call the call to set
     */
    public void setCall(String call) {
        this.call = call;
    }

    /**
     * @return the qual
     */
    public String getQual() {
        return qual;
    }

    /**
     * @param qual the qual to set
     */
    public void setQual(String qual) {
        this.qual = qual;
    }

    /**
     * @return the key
     */
    public String getKey() {
        return key;
    }

    /**
     * @param key the key to set
     */
    public void setKey(String key) {
        this.key = key;
    }

    @Override
    public String toString() {
        return this.key + ":" + this.ref + "/" + this.call + "/" + this.qual;
    }

    public final void getCharactorVector() {
        char[] result = new char[this.qual.length()];
        int index = 0;
        for (int i = 0; i < this.call.length(); i++) {
            char x = this.call.charAt(i);
            char b;
            if (x == '.' || x == ',') {
                b = '.';
            } else if (x == '+' || x == '-') {
                int j = i + 1;
                char x1 = this.call.charAt(j);
                String number = "";
                while (x1 <= '9' && x1 >= '0') {
                    number += x1;
                    j += 1;
                    x1 = this.call.charAt(j);
                }
                int inssize = Integer.valueOf(number);
                i = j + inssize - 1;
                continue;
            } else if (x == '$') {
                continue;
            } else if (x == '^') {
                i += 1;
                continue;
            } else {
                b = x;
                this.setHasb(true);
            }
            try {
                result[index] = b;
            } catch (ArrayIndexOutOfBoundsException e) {
                System.out.println(this.raw);
            }

            index += 1;
        }
        this.charactorv = result;
    }

    public final void eliminateMimorAlleles() {
        char[] chv = this.charactorv;
        this.maskv = new boolean[chv.length];
        //countv = {ref, a, t, g, c}
        int[] countv = {0, 0, 0, 0, 0};
        char[] nucvector1 = {'.', 'a', 't', 'g', 'c'};
        char[] nucvector2 = {'.', 'A', 'T', 'G', 'C'};
        int length = chv.length;
        for (int i = 0; i < length; i++) {
            if (chv[i] == '.') {
                countv[0] += 1;
            } else if (chv[i] == 'a' || chv[i] == 'A') {
                countv[1] += 1;
            } else if (chv[i] == 't' || chv[i] == 'T') {
                countv[2] += 1;
            } else if (chv[i] == 'g' || chv[i] == 'G') {
                countv[3] += 1;
            } else if (chv[i] == 'c' || chv[i] == 'C') {
                countv[4] += 1;
            }
        }
        int maxpos = 0;
        int max = -1;
        for (int i = 0; i < 5; i++) {
            if (countv[i] > max) {
                max = countv[i];
                maxpos = i;
            }
        }
        int secondpos = 0;
        int second = -1;
        for (int i = 0; i < 5; i++) {
            if (countv[i] > second && i != maxpos) {
                secondpos = i;
                second = countv[i];
            }
        }
        int truelength = 0;
        for (int i = 0; i < length; i++) {
            if (chv[i] == nucvector1[maxpos] || chv[i] == nucvector2[maxpos] || chv[i] == nucvector1[secondpos] || chv[i] == nucvector2[secondpos]) {
                maskv[i] = true;
                truelength += 1;
            } else {
                maskv[i] = false;
            }
        }
        this.setMajor(nucvector1[maxpos]);
        this.setMinor(nucvector1[secondpos]);
        this.modifiedLength = truelength;
    }

    public void printBooleanVector() {
        boolean[] temp = this.toBooleanVector();
        for (int i = 0; i < temp.length; i++) {
            if (temp[i]) {
                System.out.print("A");
            } else {
                System.out.print("B");
            }
        }
        System.out.println("");
    }

    public void printErrorVector() {
        float[] temp = this.toErrorVector();
        for (int i = 0; i < temp.length; i++) {
            System.out.print(String.valueOf(temp[i]));
            System.out.print(", ");
        }
        System.out.println("");
    }

    public double getSimpleAF() {
        boolean[] bvec = this.toBooleanVector();
        int x = bvec.length;
        int cnt = 0;
        for (int i = 0; i < x; i++) {
            if (bvec[i] == true) {
                cnt += 1;
            }
        }
        return (double) cnt / (double) x;
    }

    public double getSimpleMatchNumber() {
        boolean[] bvec = this.toBooleanVector();
        int x = bvec.length;
        int cnt = 0;
        for (int i = 0; i < x; i++) {
            if (bvec[i] == true) {
                cnt += 1;
            }
        }
        return cnt;
    }

    public boolean[] toBooleanVector() {
        boolean[] result = new boolean[this.modifiedLength];
        int orilength = this.charactorv.length;
        int index = 0;
        for (int i = 0; i < orilength; i++) {
            if (this.maskv[i] == true) {
                if (this.charactorv[i] == '.' || this.charactorv[i] == ',') {
                    result[index] = true;
                } else {
                    result[index] = false;
                }
                index += 1;
            }
        }
//        if (this.charactorv.length != this.modifiedLength) {
//            System.out.println("\n2nd minor found at " + String.valueOf(this.pos) + "\tMajor=" + this.major +", minor=" +this.minor);
//            VirmidUtil.printCharVector(this.charactorv);
//            System.out.println("");
//            VirmidUtil.printBooleanVector(this.maskv);
//            System.out.println("");
//            VirmidUtil.printBooleanVector(result);
//            System.out.println("");
//            VirmidUtil.printCharVector(this.qual.toCharArray());
//            System.out.println("");
//            VirmidUtil.printDoubleVector(this.toErrorVector());
//        }
        return result;
    }

    public float[] toErrorVector() {
        float[] error = new float[this.charactorv.length];
        for (int i = 0; i < this.qual.length(); i++) {
            char x = this.qual.charAt(i);
            float q = -(float) ((int) x - 33) / 10;
            float e = (float) Math.pow(10, q);
            error[i] = e;
        }
        float[] result = new float[this.modifiedLength];
        int orilength = this.charactorv.length;
        int index = 0;
        for (int i = 0; i < orilength; i++) {
            if (this.maskv[i] == true) {
                result[index] = error[i];
                index += 1;
            }
        }
        return result;
    }

    public float[] toMappingErrorVector() {
        float[] mappingerror = new float[this.charactorv.length];
        for (int i = 0; i < this.mappingqual.length(); i++) {
            char x = this.mappingqual.charAt(i);
            float q = -(float) ((int) x - 33) / 10;
            float e = (float) Math.pow(10, q);
            mappingerror[i] = e;
        }
        float[] result = new float[this.modifiedLength];
        int orilength = this.charactorv.length;
        int index = 0;
        for (int i = 0; i < orilength; i++) {
            if (this.maskv[i] == true) {
                result[index] = mappingerror[i];
                index += 1;
            }
        }
        return result;
    }

    public double getRegionMappingQulity() {
        float[] mappingerror = this.toMappingErrorVector();
        double ave = 0.0;
        for (int i = 0; i < mappingerror.length; i++) {
            ave += mappingerror[i];
        }
        ave = ave / mappingerror.length;
        return ave;
    }

    public double getRegionMaximumMappingQuality() {
        float[] mappingerror = this.toMappingErrorVector();
        double ave = 0.0;
        for (int i = 0; i < mappingerror.length; i++) {
            if (ave < mappingerror[i]) {
                ave = mappingerror[i];
            }
        }
        return ave;
    }

    public double getRegionMeanMappingQuality() {
        float[] mappingerror = this.toMappingErrorVector();
        double ave = 0.0;
        for (int i = 0; i < mappingerror.length; i++) {
            ave += mappingerror[i];
        }
        return ave / mappingerror.length;
    }

    public int getRegionMeanMappingPhred() {
        int totalphred = 0;
        for (int i = 0; i < this.mappingqual.length(); i++) {
            char x = this.mappingqual.charAt(i);
            totalphred += (int) x - 33;
        }
        return totalphred / mappingqual.length();
    }

    /**
     * @return the valid
     */
    public boolean isValid() {
        return valid;
    }

    /**
     * @param valid the valid to set
     */
    public void setValid(boolean valid) {
        this.valid = valid;
    }

    /**
     * @return the major
     */
    public char getMajor() {
        return major;
    }

    /**
     * @param major the major to set
     */
    public void setMajor(char major) {
        this.major = major;
    }

    /**
     * @return the minor
     */
    public char getMinor() {
        return minor;
    }

    /**
     * @param minor the minor to set
     */
    public void setMinor(char minor) {
        this.minor = minor;
    }

    /**
     * @return the hasb
     */
    public boolean hasB() {
        return hasb;
    }

    /**
     * @param hasb the hasb to set
     */
    public void setHasb(boolean hasb) {
        this.hasb = hasb;
    }
}
