package hybridsearch;

import java.util.HashMap;

public class ConditionalProbability {
	static private HashMap<ConditionalProbability, Float> conditionalProbMap = new HashMap<ConditionalProbability, Float>(); 
	static private float minVal = 1e-12f;  
	//static private int libMaxScore;
	//static private int libMinScore;
	//static private int dbMaxScore;
	//static private int dbMinScore;
	
	private int libScore;
	private int dbScore;
	private float prob;
	
	protected ConditionalProbability(int dbScore, int libScore){ // for training
		this.dbScore = dbScore;	
		this.libScore = libScore;
		this.prob = retrieve();
	}
	
	protected ConditionalProbability(String s){
		String[] token = s.split("\t");
		
		this.dbScore = Integer.parseInt(token[0]);	
		this.libScore = Integer.parseInt(token[1]);
		this.prob = Float.parseFloat(token[2]);
	}
	
	public int getDBScore(){
		return dbScore;
	}
	
	public int getLibScore(){
		return libScore;
	}
	
	private float retrieve(){
		Float v = conditionalProbMap.get(this);
		if(v == null) return minVal;// TODO
		else return v;
	}
	
	public float get(){
		return prob;
	}
	
	protected void set(float p){
		conditionalProbMap.put(this, p);
	}
	
	public String toString(){
		return dbScore + "\t" + libScore + "\t" + get();
	}
	
	public int hashCode(){
		return libScore +  dbScore << 10;
	}
	
	public boolean equals(Object other){
		if(other instanceof ConditionalProbability){
			ConditionalProbability o = (ConditionalProbability)other;
			
			return o.dbScore == this.dbScore && o.libScore == this.libScore;
		}
		
		return false;
	}
	
	
	
	
}
