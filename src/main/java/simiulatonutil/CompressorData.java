package simiulatonutil;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.27 027 14:40:54
 */
public class CompressorData {
    public static void main(String[] args) throws IOException {

        List<String> data = FileUtils.readLines(new File("CompressorData.txt"), "UTF-8");
        int len = data.size() - 1;
        String[] r = data.get(0).split("\\s");
        ArrayList<Double> head_m1 = new ArrayList<>();
        ArrayList<Double> head_m2 = new ArrayList<>();
        ArrayList<Double> efficiency1 = new ArrayList<>();
        ArrayList<Double> efficiency2 = new ArrayList<>();
        ArrayList<Double> flowrate1 = new ArrayList<>();
        ArrayList<Double> flowrate2 = new ArrayList<>();
        ArrayList<Double> pressratio1 = new ArrayList<>();
        ArrayList<Double> pressratio2 = new ArrayList<>();
        ArrayList<Double> rotationrate1 = new ArrayList<>();
        ArrayList<Double> rotationrate2 = new ArrayList<>();
        data.stream().skip(1).forEach(d -> {
            String[] line = d.split("\\s");
            head_m1.add(Double.parseDouble(line[0]));
            efficiency1.add(Double.parseDouble(line[1]));
            flowrate1.add(Double.parseDouble(line[2]));
            pressratio1.add(Double.parseDouble(line[3]));
            rotationrate1.add(Double.parseDouble(line[4]));
            head_m2.add(Double.parseDouble(line[5]));
            efficiency2.add(Double.parseDouble(line[6]));
            flowrate2.add(Double.parseDouble(line[7]));
            pressratio2.add(Double.parseDouble(line[8]));
            rotationrate2.add(Double.parseDouble(line[9]));
        });
        head_m1.addAll(head_m2);
        efficiency1.addAll(efficiency2);
        flowrate1.addAll(flowrate2);
        pressratio1.addAll(pressratio2);
        rotationrate1.addAll(rotationrate2);

        double[] head_m = head_m1.stream().mapToDouble(Double::doubleValue).toArray();
        double[] efficiency = efficiency1.stream().mapToDouble(Double::doubleValue).toArray();
        double[] flowrate_am3_h = flowrate1.stream().mapToDouble(Double::doubleValue).toArray();
        double[] pressratio = pressratio1.stream().mapToDouble(Double::doubleValue).toArray();
        double[] rotationrate_RPM = rotationrate1.stream().mapToDouble(Double::doubleValue).toArray();

        WeightedObservedPoints obs = new WeightedObservedPoints();
        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(1);

        double test = rotationrate_RPM[0];
        for (int i = 0; i < rotationrate_RPM.length; i++) {
            if (rotationrate_RPM[i] != test) {
                double[] coeff = fitter.fit(obs.toList());
                PolynomialFunction f = new PolynomialFunction(coeff);
                System.out.println("RPM = " + test + "\t" + f);
                obs.clear();
                test = rotationrate_RPM[i];
            }
            double x = flowrate_am3_h[i];
            double y = pressratio[i];
            obs.add(x * x, y * y);
        }
        double[] coeff = fitter.fit(obs.toList());
        PolynomialFunction f = new PolynomialFunction(coeff);
        System.out.println("RPM = " + test + "\t" + f);
        obs.clear();

        test = rotationrate_RPM[0];
        System.out.print("min = " + flowrate_am3_h[0]);
        for (int i = 0; i < rotationrate_RPM.length; i++) {
            if (rotationrate_RPM[i] != test) {
                System.out.println("\tmax = " + flowrate_am3_h[i - 1]);
                System.out.print("min = " + flowrate_am3_h[i]);
                test = rotationrate_RPM[i];
            }
        }
        System.out.println("\tmax = " + flowrate_am3_h[flowrate_am3_h.length - 1]);
    }
}
