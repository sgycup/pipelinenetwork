package simiulatonutil;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.19 019 11:29:29
 */
public class CapacityCalculator {
    BWRS bwrs;
    double B, C, D, E, F, R;//R单位为kJ/(kg*K)
    public CapacityCalculator(double[] parameters, BWRS bwrs) {
        B = parameters[0];
        C = parameters[1];
        D = parameters[2];
        E = parameters[3];
        F = parameters[4];
        R = parameters[5];
        this.bwrs = bwrs;
    }

    /**
     * @Author Shi Guoyun
     * @Description 计算定压比容
     * @Date 11:47:16 2021.01.19 019
     *
     * @param density 密度，kmol/m3
     * @param temperature 温度，K
     * @return kJ/(kg*K)
     */
    public double capacityConstantPressure(double density, double temperature) {
        double cv0 = B + 2 * C * temperature + 3 * D * Math.pow(temperature, 2) + 4 * E * Math.pow(temperature, 3)
                + 5 * F * Math.pow(temperature, 4) - R;
        double Mg = bwrs.parameter.Mg;
        double deltacv = bwrs.deltaCapacity(density, temperature) /Mg;
        double cv = cv0 + deltacv;
//        double para1 = Math.pow(bwrs.partialP_Temperature(density, temperature), 2);
//        double para2 = bwrs.partialP_Density(density, temperature);
//        double cp = cv + temperature / Math.pow(density, 2) * para1 / para2 / Mg;
        return cv;
    }
}
