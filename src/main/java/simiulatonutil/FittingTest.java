package simiulatonutil;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.26 026 14:36:06
 */
public class FittingTest {
    public static void main(String[] args) {

        final WeightedObservedPoints obs = new WeightedObservedPoints();
        obs.add(2,  0.7 * 2 + 20 + 0.4);
        obs.add(12,  0.7 * 12 + 20 + 0.3);
        obs.add(32,  0.7 * 32 + 20 + 3.4);
        obs.add(34 ,  0.7 * 34 + 20 + 5.8);
        obs.add(58 , 0.7 * 58 + 20 + 8.4);
        obs.add(43 , 0.7 * 43 + 20 + 0.28);
        obs.add(27 , 0.7 * 27 + 20 + 0.4);

        // Instantiate a two-degree polynomial fitter.
        final PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2);

        // Retrieve fitted parameters (coefficients of the polynomial function).
        final double[] coeff = fitter.fit(obs.toList());
        for (double c : coeff) {
            System.out.println(c);
        }
        PolynomialFunction f = new PolynomialFunction(coeff);
        System.out.println(f);
    }
}
