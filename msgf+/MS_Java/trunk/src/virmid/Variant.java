/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package virmid;

/**
 *
 * @author rachmani
 */
public class Variant {
    private String chrom;
    private int pos;
    private String key;
    private String id;
    private String ref;
    private String alt;
    private float qual;
    private String filter;
    private String info;
    private String format;
    private String[] annot;
    private String raw;
    private boolean valid = false;
    
    public Variant() {
        init();
    }
    
    public Variant(String raw) {
        init();
        String[] fields = raw.split("\t");
        if (fields.length >= 7) {
            this.chrom = fields[0];
            this.pos = Integer.valueOf(fields[1]);
            this.id = fields[2];
            this.ref = fields[3];
            this.alt = fields[4];
            this.qual = Float.valueOf(fields[5]);
            this.filter = fields[6];
            this.info = fields[7];
            this.format = fields[8];
            this.key = this.chrom + ":" + String.valueOf(this.pos);
            this.raw = raw;
            this.valid = true;
        }
    }

    private void init() {
        this.setChrom("");
        this.setPos(0);
        this.setId("");
        this.setRef("");
        this.setAlt("");
        this.setQual(0);
        this.setFilter("");
        this.setInfo("");
        this.setInfo("");
        this.setFormat("");
        
                
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
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(String id) {
        this.id = id;
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
     * @return the alt
     */
    public String getAlt() {
        return alt;
    }

    /**
     * @param alt the alt to set
     */
    public void setAlt(String alt) {
        this.alt = alt;
    }

    /**
     * @return the qual
     */
    public float getQual() {
        return qual;
    }

    /**
     * @param qual the qual to set
     */
    public void setQual(float qual) {
        this.qual = qual;
    }

    /**
     * @return the filter
     */
    public String getFilter() {
        return filter;
    }

    /**
     * @param filter the filter to set
     */
    public void setFilter(String filter) {
        this.filter = filter;
    }

    /**
     * @return the info
     */
    public String getInfo() {
        return info;
    }

    /**
     * @param info the info to set
     */
    public void setInfo(String info) {
        this.info = info;
    }

    /**
     * @return the format
     */
    public String getFormat() {
        return format;
    }

    /**
     * @param format the format to set
     */
    public void setFormat(String format) {
        this.format = format;
    }

    /**
     * @return the annot
     */
    public String[] getAnnot() {
        return annot;
    }

    /**
     * @param annot the annot to set
     */
    public void setAnnot(String[] annot) {
        this.annot = annot;
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
        return this.chrom + ":" + String.valueOf(this.pos) + ":" + this.ref + "->" + this.alt + " (" + 
                String.valueOf(this.qual) + ")";
    }
    
    public String toBED() {
        return this.chrom + "\t" + String.valueOf(pos-1) + "\t" + String.valueOf(pos);
    }
    
    public String toVCFRecord() {
        return this.raw;
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
}
