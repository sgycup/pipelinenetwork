package simiulatonutil;

import cnic.pipesimulation.SinglePipeDiscretizer;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.13 013 9:40:11
 */
public class PropertyUpdater {
    public SinglePipeDiscretizer discretizer;
    public BWRS bwrs;
    public CapacityCalculator capacity;
    public Colebrook colebrook;

    public void propertyUpdate() {
        int nump = discretizer.nump;
        double[] c2 = discretizer.c2;
        double[] density = discretizer.rho;
        double[] pressure = discretizer.p;
        double[] temperature = discretizer.T;
        double[] cp = discretizer.Cp;
        double Mg = bwrs.parameter.Mg;
//        // 更新密度
//        for (int i = 0; i < nump; i++) {
//            density[i] = bwrs.getDensity(pressure[i] / 1000, temperature[i]) * Mg;
//        }
        // 更新声速
        for (int i = 0; i < nump; i++) {
            double t = temperature[i];
            c2[i] = bwrs.partialP_Density(density[i] / Mg, t) * 1000 / Mg;
        }
        // 更新热容
        for (int i = 0; i < nump; i++) {
            cp[i] = capacity.capacityConstantPressure(density[i] / Mg, temperature[i]) * 1000;
        }
        int numv = discretizer.numv;
        double D = discretizer.D;
        double[] u = discretizer.u;
        double[] lambda = discretizer.lambda;
        // 更新摩阻系数
        for (int i = 0; i < numv; i++) {
            double rhoB = (density[i] + density[i + 1]) / 2;
            double tB = (temperature[i] + temperature[i + 1]) / 2;
            double mu = ParameterCalculator.viscosity(rhoB, tB);
            double Re = rhoB * Math.abs(u[i]) * D / mu;
            if (Double.isNaN(Re)) {
                throw new RuntimeException("Re is NaN (rho = " + rhoB + ", u = " + u[i] + ", mu = " + mu + ")");
            }
            lambda[i] = ParameterCalculator.lambdaColebrook(colebrook, Re);
        }
    }
}
