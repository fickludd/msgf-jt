package IMS;

import java.util.ArrayList;
import java.util.HashMap;

import msutil.IonType;

public class IonPairCounter {
	private static HashMap<FragmentParameter, HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>> ionCountMapTarget = new HashMap<FragmentParameter, HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>>();
	private static  HashMap<FragmentParameter, HashMap<IonType, Float>> totalTargetCount =  new HashMap<FragmentParameter, HashMap<IonType, Float>>();
	private static HashMap<FragmentParameter, HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>> ionCountMapDecoy = new HashMap<FragmentParameter, HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>>();
	private static  HashMap<FragmentParameter, HashMap<IonType, Float>> totalDecoyCount =  new HashMap<FragmentParameter, HashMap<IonType, Float>>();
	private static boolean isTargetNormalized = false, isDecoyNormalized = false; 
	
	public static float getProbability(FragmentParameter para, IonType ion1, IonType ion2, int ratio, boolean isTarget){
		HashMap<FragmentParameter, HashMap<IonType, Float>> total;
		if(ion1.equals(ion2) || ion1.getCharge()!=ion2.getCharge()) return 0;		
		
		HashMap<FragmentParameter, HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>> map;
		if(isTarget){
			if(!isTargetNormalized) normalize(true);
			map = ionCountMapTarget;
			total = totalTargetCount;
		}
		else{
			if(!isDecoyNormalized) normalize(false);
			map = ionCountMapDecoy;
			total = totalDecoyCount;
		}
		
		if(map.containsKey(para)){
			HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>> ionCounter = map.get(para);
			if(ionCounter.containsKey(ion1)){
				HashMap<IonType, HashMap<Integer, Float>> subIonCounter = ionCounter.get(ion1);
				if(subIonCounter.containsKey(ion2)){
					HashMap<Integer, Float> subsubIonCounter = subIonCounter.get(ion2);
					if(subsubIonCounter.containsKey(ratio)) return subsubIonCounter.get(ratio);
					else return .5f/total.get(para).get(ion1);
				}else return 0;
			}else return 0;
		}else return 0;		
	}
	
	private static void normalize(boolean isTarget){
		HashMap<FragmentParameter, HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>> map;
		HashMap<FragmentParameter, HashMap<IonType, Float>> total;
		if(isTarget){
			map = ionCountMapTarget;
			total = totalTargetCount;
			isTargetNormalized = true;
		}
		else{
			map = ionCountMapDecoy;
			total = totalDecoyCount;
			isDecoyNormalized = true;
		}
		
		for(FragmentParameter para : map.keySet()){
			HashMap<IonType, Float> t = total.get(para);			
			HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>> ionCounter = map.get(para);
			
			for(IonType ion1 : ionCounter.keySet()){
				float it = t.get(ion1);
				HashMap<IonType, HashMap<Integer, Float>> subIonCounter = ionCounter.get(ion1);
				for(IonType ion2: subIonCounter.keySet()){
					if(ion1.equals(ion2) || ion1.getCharge()!=ion2.getCharge()) continue;
					HashMap<Integer, Float> subsubIonCounter = subIonCounter.get(ion2);
					for(int k : subsubIonCounter.keySet()){
						subsubIonCounter.put(k, (subsubIonCounter.get(k)+1)/(it+1));
					}
				}
			}			
		}
	}
	
	public static void update(FragmentParameter para, HashMap<IonType, Float> ions, boolean isTarget){
		HashMap<FragmentParameter, HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>> map;
		HashMap<FragmentParameter, HashMap<IonType, Float>> total;
		if(isTarget){
			map = ionCountMapTarget;
			total = totalTargetCount;
		}
		else{
			map = ionCountMapDecoy;
			total = totalDecoyCount;
		}
		
		if(!map.containsKey(para)){
			map.put(para, new HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>>());
			total.put(para, new HashMap<IonType, Float>());
		}
		
		HashMap<IonType, HashMap<IonType, HashMap<Integer, Float>>> ionCounter = map.get(para);
		HashMap<IonType, Float> ionTotal = total.get(para);
		
		for(IonType ion1 : ions.keySet()){
			if(!ionCounter.containsKey(ion1)){
				ionCounter.put(ion1, new HashMap<IonType, HashMap<Integer, Float>>());
				ionTotal.put(ion1, 0f);				
			}
			HashMap<IonType, HashMap<Integer, Float>> subIonCounter = ionCounter.get(ion1);
			for(IonType ion2 : ions.keySet()){
				if(ion1.equals(ion2) || ion1.getCharge()!=ion2.getCharge()) continue;
				if(!subIonCounter.containsKey(ion2)){
					subIonCounter.put(ion2, new HashMap<Integer, Float>());
				}
				HashMap<Integer, Float> subsubIonCounter = subIonCounter.get(ion2);
				int ratio = getRatio(ions.get(ion1), ions.get(ion2));
				if(!subsubIonCounter.containsKey (ratio)){
					subsubIonCounter.put(ratio, 0f);
				}
				subsubIonCounter.put(ratio, subsubIonCounter.get(ratio)+1);
								
			}
			ionTotal.put(ion1, ionTotal.get(ion1)+1);			
		}		
	}
	
	private static int getRatio(float v1, float v2){
		if(v1 == 0){
			if(v2 == 0) return -11;
			else return -12;
		}else{
			if(v2 == 0) return 11;
			else{
				float r = 0, f = 1;;
				if(v1 > v2){
					r = v1/v2;
				}else{
					r = v2/v1;
					f=-1;
				}				
				return (int)(Math.min(10, r) * f);
			}
		}
	}
	
	public static ArrayList<Integer> getAllRatios(){
		ArrayList<Integer> ret = new ArrayList<Integer>();
		for(int i=-12;i<12;i++)
			ret.add(i);
		return ret;
	}
	
	public static void main(String[] args){
		for(int v1=0;v1<20;v1++){
			for(int v2=0;v2<20;v2++){
				System.out.println(v1+" " + v2 + " " + getRatio(v1,v2));
			}	
		}
	}
	
}
