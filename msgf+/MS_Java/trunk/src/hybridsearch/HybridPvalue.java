package hybridsearch;


public class HybridPvalue {
	static private float alpha=1;
	static private float delta = 0.25f;
	static final public int maxScore = (int)-Math.log10(1e-17f); // TODO maxScore should guarantee 100 PSMs... 
	//static private HashMap<ArrayList<Float>, Float> cpTable = new HashMap<ArrayList<Float>, Float>();
	
	private float specProb, slgfPvalue=-Float.NaN;
	 
	public HybridPvalue(float specProb, float slgfPvalue) {
		this.specProb = Math.max(1e-40f, specProb);
		this.slgfPvalue = slgfPvalue;
	}
	
	public HybridPvalue(float specProb) {
		this.specProb = Math.max(1e-40f, specProb);
	}
	
	public float get(){
		float p;
		if(Float.isNaN(slgfPvalue)) 
			p = getExpectedPvalue(getDBScore());
		else p = getAccuratePvalue(getDBScore(), getLibScore());
		
		//System.out.println(specProb + "\t" + slgfPvalue + "\t" + p);
		return p;
	}
	
	private float getDBScore(){
		return (float) -Math.log10(specProb);
	}
	
	private float getLibScore(){
		float z = 4;
		float t = (float) (1f/(slgfPvalue + Math.pow(10, -z)) - slgfPvalue/(1+2*Math.pow(10, -z)));
		float out = (float) ((-Math.log10(t) + z) * maxScore/2/z  );
			
			//(float) Math.log10(slgfPvalue*1e3f)*maxScore/3;
		//if(Float.isNaN(out))
		//	System.out.println(slgfPvalue + "\t" + out);
		return out;
		//return (float) slgfPvalue * maxScore;// -(Math.log10(1-slgfPvalue)/Math.log10(1.0665));//TODO adjust
	}
	
	protected int getIntegerDBScore(){
		return (int) Math.min(maxScore-1, Math.round(getDBScore()));
	}
	
	protected int getIntegerLibScore(){
		return (int) Math.min(maxScore-1, Math.round(getLibScore()));//TODO is this good? .. or..
	}
	
	private float getConditionalProbability(float b, float a){
	//	ArrayList<Float> key = new ArrayList<Float>();
		int br = Math.round(b);
		
	//	key.add((float)br); key.add(a);
		
	//	if(cpTable.containsKey(key)){
	//		return cpTable.get(key);
	//	}
		
		float p = 0;
		
		for(float k=a;k<maxScore;k++){
			int kf = (int)k;
			int kc = kf + 1;
			ConditionalProbability cp1 = new ConditionalProbability(br, kf);
			ConditionalProbability cp2 = new ConditionalProbability(br, kc);
			p += Math.max(0, (kc-k) * cp1.get() + (k-kf) * cp2.get());
			//System.out.println(cp1.get() + "\t" + cp2.get() + "\t" + kf + "\t" + kc);
			
		}
		
	//	cpTable.put(key, p);
		return p;
	}
	//TODO how to manage maxScore??
	
	private float getAccuratePvalue(float db, float lib){
		float t = db + alpha * lib;
		float p = (float) Math.pow(10, -t);
		
		for(int k=0; k<Math.floor(t/delta-1);k++){
			float f = (float) (Math.pow(10, -k*delta) - Math.pow(10, -(k+1)*delta));
			p += f * getConditionalProbability(k*delta, (t-k*delta)/alpha);
		}
		return p;
	}
	
	private float getExpectedPvalue(float db){
		float p = 0;
		for(int l=0; l<(maxScore-1)/delta; l++){
			float x = getConditionalProbability(db, delta*l) - getConditionalProbability(db, delta*(l+1));
			p += getAccuratePvalue(db, delta * l) * Math.max(0, x);
		}
		return p;
	}
	
	static protected void setAlpha(float a){
		alpha = a;
	}
	
}
