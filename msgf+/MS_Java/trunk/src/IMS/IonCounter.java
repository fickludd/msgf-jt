package IMS;

import java.util.HashMap;

import msutil.IonType;

public class IonCounter {
	private static HashMap<FragmentParameter, HashMap<IonType, Float>> ionCountMapTarget = new HashMap<FragmentParameter, HashMap<IonType, Float>>();
	private static  HashMap<FragmentParameter, Integer> totalTargetCount =  new HashMap<FragmentParameter, Integer>();
	private static HashMap<FragmentParameter, HashMap<IonType, Float>> ionCountMapDecoy = new HashMap<FragmentParameter, HashMap<IonType, Float>>();
	private static HashMap<FragmentParameter, Integer> totalDecoyCount =  new HashMap<FragmentParameter, Integer>();
	private static boolean isTargetNormalized = false, isDecoyNormalized = false; 
	
	public static float getProbability(FragmentParameter para, IonType ion, boolean isTarget){
		HashMap<FragmentParameter, HashMap<IonType, Float>> map;
		HashMap<FragmentParameter, Integer> total;
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
			if(map.get(para).containsKey(ion)){
				return map.get(para).get(ion);
			}else return .5f/total.get(para);
		}else return .5f/total.get(para);
	}
	
	private static void normalize(boolean isTarget){
		HashMap<FragmentParameter, HashMap<IonType, Float>> map;
		HashMap<FragmentParameter, Integer> total;
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
			int t = total.get(para);
			HashMap<IonType, Float> ionCounter = map.get(para);
			for(IonType ion : ionCounter.keySet()){
				ionCounter.put(ion, (ionCounter.get(ion)+1)/(t+1));
			}
		}
	}
	
	public static void update(FragmentParameter para, HashMap<IonType, Float> ions, boolean isTarget){
		HashMap<FragmentParameter, HashMap<IonType, Float>> map;
		HashMap<FragmentParameter, Integer> total;
		if(isTarget){
			map = ionCountMapTarget;
			total = totalTargetCount;
		}
		else{
			map = ionCountMapDecoy;
			total = totalDecoyCount;
		}
		
		if(!map.containsKey(para)){
			map.put(para, new HashMap<IonType, Float>());
			total.put(para, 0);
		}
		HashMap<IonType, Float> ionCounter = map.get(para);		
		for(IonType ion : ions.keySet()){
			if(ions.get(ion) <=0) continue;
			Float count = ionCounter.get(ion);
			if(count == null) count = 0f;			
			count ++;
			ionCounter.put(ion, count);
		}
		total.put(para, total.get(para)+1);
	}
}
