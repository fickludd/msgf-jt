package uninovoOld.features;

import java.util.Comparator;

import msgf.Tolerance;

public class FeatureComparator implements Comparator<Feature>{ // 
	private Tolerance tol;
	public FeatureComparator(Tolerance tol){
		this.tol = tol;
	}
	@Override
	public int compare(Feature arg0, Feature arg1) {
		return 0;//Float.compare(arg0.getH(tol), arg1.getH(tol));
		
		//return Float.compare(Feature.getKLDivergenceFromNullCondition(ion, arg0), Feature.getKLDivergenceFromNullCondition(ion, arg1));
	}

}
