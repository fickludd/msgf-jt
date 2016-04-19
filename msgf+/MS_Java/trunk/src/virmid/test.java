package virmid;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;

public class test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		double a = 1.1234;
		double b = 1e-1;
		
		System.out.println(Math.log(a));

		System.out.println(Math.log(b));
		
		System.out.println(new BigDecimal(a).add(new BigDecimal(b)).equals(new BigDecimal(a)));
		
		long t = System.nanoTime();
		double c = 0;
		for(float i=0; i<10000000;i++){
			c *= a;
		}
		System.out.println(new BigDecimal("1e-50000").doubleValue() + "\t" + 
				
				"done" + c + "\n" + (System.nanoTime()-t));
		
		 t = System.nanoTime();

		BigDecimal bc = new BigDecimal(1, new MathContext(100));
		BigDecimal ba = new BigDecimal(a);
		for(float i=0; i<100;i++){
			bc = bc.multiply(ba,  new MathContext(100));
		}
		ba.abs();
		
	
		System.out.println("done" + bc.precision() + "\n" + (System.nanoTime()-t));
		System.out.println((new BigDecimal(100.8).scale()));
	}

	public static BigDecimal log(BigDecimal x) {
        BigDecimal result = BigDecimal.ZERO;

        BigDecimal input = new BigDecimal(x.toString());
        int decimalPlaces = 500;
        int scale = input.precision() + decimalPlaces;

        int maxite = 10000;
        int ite = 0;
        BigDecimal maxError_BigDecimal = new BigDecimal(BigInteger.ONE,decimalPlaces + 1);
        //System.out.println("maxError_BigDecimal " + maxError_BigDecimal);
       // System.out.println("scale " + scale);

        RoundingMode a_RoundingMode = RoundingMode.UP;

        BigDecimal two_BigDecimal = new BigDecimal("2");
        BigDecimal base_BigDecimal = new BigDecimal(Math.E);

        while (input.compareTo(base_BigDecimal) == 1) {
            result = result.add(BigDecimal.ONE);
            input = input.divide(base_BigDecimal, scale, a_RoundingMode);
        }

        BigDecimal fraction = new BigDecimal("0.5");
        input = input.multiply(input);
        BigDecimal resultplusfraction = result.add(fraction);
        while (((resultplusfraction).compareTo(result) == 1)
                && (input.compareTo(BigDecimal.ONE) == 1)) {
            if (input.compareTo(base_BigDecimal) == 1) {
                input = input.divide(base_BigDecimal, scale, a_RoundingMode);
                result = result.add(fraction);
            }
            input = input.multiply(input);
            fraction = fraction.divide(two_BigDecimal, scale, a_RoundingMode);
            resultplusfraction = result.add(fraction);
            if (fraction.abs().compareTo(maxError_BigDecimal) == -1){
                break;
            }
            if (maxite == ite){
                break;
            }
            ite ++;
        }

        MathContext a_MathContext = new MathContext(((decimalPlaces - 1) + (result.precision() - result.scale())),RoundingMode.HALF_UP);
        BigDecimal roundedResult = result.round(a_MathContext);
        BigDecimal strippedRoundedResult = roundedResult.stripTrailingZeros();
        //return result;
        //return result.round(a_MathContext);
        return strippedRoundedResult;
    }
	
}
