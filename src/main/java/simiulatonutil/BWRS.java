package simiulatonutil;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.12 012 13:48:20
 */
public class BWRS {
    public final static double[] Ai = {0.44369, 1.28438, 0.356306, 0.544979, 0.528629, 0.484011, 0.070523, 0.504087, 0.030745, 0.073283, 0.00645};
    public final static double[] Bi = {0.115449, -0.920731, 1.70871, -0.270896, 0.349261, 0.75413, -0.044448, 1.32245, 0.179433, 0.463492, -0.02214};

    public final ComponentParameter parameter;
    public final double R = 8.3143;

    public BWRS(ComponentParameter parameter) {
        this.parameter = parameter;
    }

    /**
     * @param density 密度，kmol/m3
     * @param temperature 温度，K
     * @return 压力，kPa
     */
    public double value(double density, double temperature) {
        double p = density * R * temperature;
        p += (parameter.B0 * R * temperature - parameter.A0 - parameter.C0 / Math.pow(temperature, 2) + parameter.D0
                / Math.pow(temperature, 3) - parameter.E0 / Math.pow(temperature, 4)) * Math.pow(density, 2);
        p += (parameter.b * R * temperature - parameter.a - parameter.d / temperature) * Math.pow(density, 3);
        p += parameter.alpha * (parameter.a + parameter.d / temperature) * Math.pow(density, 6);
        p += parameter.c * Math.pow(density, 3) / Math.pow(temperature, 2) * (1 + parameter.gamma * Math.pow(density, 2))
                * Math.exp(-parameter.gamma * Math.pow(density, 2));
        return p;
    }

    /**
     * @Author Shi Guoyun
     * @Description 压力对密度的偏导数
     * @Date 15:00:19 2021.01.12 012
     *
     * @param density 密度，kmol/m3
     * @param temperature 温度，K
     * @return kPa/(kmol/m3)
     */
    public double partialP_Density(double density, double temperature) {
        double dpdro = R * temperature;
        dpdro += 2 * (parameter.B0 * R * temperature - parameter.A0 - parameter.C0 / Math.pow(temperature, 2)
                + parameter.D0 / Math.pow(temperature, 3) - parameter.E0 / Math.pow(temperature, 4)) *density;
        dpdro += 3 * (parameter.b * R * temperature - parameter.a - parameter.d / temperature) * Math.pow(density, 2);
        dpdro += 6 * parameter.alpha * (parameter.a + parameter.d / temperature) * Math.pow(density, 5);
        dpdro += 3 * parameter.c * Math.pow(density, 2) / Math.pow(temperature, 2)
                * (1 + parameter.gamma * Math.pow(density, 2) - 2.0 / 3 * Math.pow(parameter.gamma, 2)
                * Math.pow(density, 4)) * Math.exp(-parameter.gamma * Math.pow(density, 2));
        return dpdro;
    }

    /**
     * @Author Shi Guoyun
     * @Description 压力对温度的偏导数
     * @Date 9:18:08 2021.01.18 018
     *
     * @param density 密度，kmol/m3
     * @param temperature 温度，K
     * @return kPa/K
     */
    public double partialP_Temperature(double density, double temperature) {
        double dpdt = density * R;
        dpdt += (parameter.B0 * R + 2 * parameter.C0 / Math.pow(temperature, 3) - 3 * parameter.D0 / Math.pow(temperature,4)
                + 4 * parameter.E0 / Math.pow(temperature, 5)) * Math.pow(density, 2);
        dpdt += (parameter.b * R + parameter.d / Math.pow(temperature, 2)) * Math.pow(density, 3);
        dpdt -= parameter.alpha * parameter.d * Math.pow(density, 6) / Math.pow(temperature, 2);
        dpdt -= 2 * parameter.c * Math.pow(density, 3) / Math.pow(temperature, 3) * (1 + parameter.gamma
                * Math.pow(density, 2)) * Math.exp(-parameter.gamma * Math.pow(density, 2));
        return dpdt;
    }

    /**
     * @Author Shi Guoyun
     * @Description 定容比热容相较于理想气体低压下高压时的差值，即delta capacity under constant volume
     * @Date 10:34:31 2021.01.19 019
     *
     * @param density 密度，kmol/m3
     * @param temperature 温度，K
     * @return kJ/(kg*K)
     */
    public double deltaCapacity(double density, double temperature) {
        double dcv = (6 * parameter.C0 / Math.pow(temperature, 3) - 12 * parameter.D0 / Math.pow(temperature, 4) + 20 * parameter.E0 / Math.pow(temperature, 5)) * density;
        dcv += parameter.d / Math.pow(temperature, 2) * Math.pow(density, 2);
        dcv += 3 * parameter.c / parameter.gamma / Math.pow(temperature, 3) * ((2 + parameter.gamma * Math.pow(density, 2)) * Math.exp(-parameter.gamma * Math.pow(density, 2)) - 2);
        return dcv;
    }

    /**
     * @Author Shi Guoyun
     * @Description 抛物线法求解密度
     * @Date 14:59:17 2021.01.12 012
     *
     * @param pressure 压力，kPa
     * @param temperature 温度，K
     * @return 密度，kmol/m3
     */
    public double getDensity(double pressure, double temperature) {
        double rho0 = 0, rho1 = pressure / R / temperature, rho2;
        double f0, f1;
        while (true) {
            f0 = value(rho0, temperature) - pressure;
            f1 = value(rho1, temperature) - pressure;
            rho2 = rho0 * f1 - rho1 * f0;
            rho2 /= f1 - f0;
            if (Math.abs(rho2 - rho1) < 1e-6) {
                rho1 = rho2;
                break;
            }
            rho0 = rho1;
            rho1 = rho2;
        }
        return  rho1;
    }

    public static class ComponentParameter {
        public double A0, B0, C0, D0, E0, a, b, c, d, alpha, gamma, Mg;
    }

    /**
     * @param criticalTemperature 临界温度，K
     * @param criticalDensity 临界密度，kmol/m3
     * @param acentricFactor 偏心因子
     * @param moleweight 摩尔质量，g/mol
     * @return BWRS方程中的11个参数
     */
    public static ComponentParameter getParameter(double criticalTemperature, double criticalDensity, double acentricFactor, double moleweight) {
        ComponentParameter parameter = new ComponentParameter();
        double R = 8.3143;
        parameter.A0 = (Ai[1] + Bi[1] * acentricFactor) * R * criticalTemperature / criticalDensity;
        parameter.B0 = (Ai[0] + Bi[0] * acentricFactor) / criticalDensity;
        parameter.C0 = (Ai[2] + Bi[2] * acentricFactor) * R * Math.pow(criticalTemperature, 3) / criticalDensity;
        parameter.D0 = (Ai[8] + Bi[8] * acentricFactor) * R * Math.pow(criticalTemperature, 3) / criticalDensity;
        parameter.E0 = (Ai[10] + Bi[10] * acentricFactor * Math.exp(-3.8 * acentricFactor)) * R * Math.pow(criticalTemperature, 5) / criticalDensity;
        parameter.a = (Ai[5] + Bi[5] * acentricFactor) * R * criticalTemperature / Math.pow(criticalDensity, 2);
        parameter.b = (Ai[4] + Bi[4] * acentricFactor) / Math.pow(criticalDensity, 2);
        parameter.c = (Ai[7] + Bi[7] * acentricFactor) * R * Math.pow(criticalTemperature, 3) / Math.pow(criticalDensity, 2);
        parameter.d = (Ai[9] + Bi[9] * acentricFactor) * R * Math.pow(criticalTemperature, 2) / Math.pow(criticalDensity, 2);
        parameter.alpha = (Ai[6] + Bi[6] * acentricFactor) / Math.pow(criticalDensity, 3);
        parameter.gamma = (Ai[3] + Bi[3] * acentricFactor) / Math.pow(criticalDensity, 2);
        parameter.Mg = moleweight;
        return  parameter;
    }

    public static ComponentParameter getMixtureParameter(ComponentParameter[] componentParameters, double[] moleFractions) {
        ComponentParameter parameter = new ComponentParameter();
        int n = moleFractions.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                parameter.A0 += moleFractions[i] * moleFractions[j] * Math.sqrt(componentParameters[i].A0 * componentParameters[j].A0);
                parameter.C0 += moleFractions[i] * moleFractions[j] * Math.sqrt(componentParameters[i].C0 * componentParameters[j].C0);
                parameter.D0 += moleFractions[i] * moleFractions[j] * Math.sqrt(componentParameters[i].D0 * componentParameters[j].D0);
                parameter.E0 += moleFractions[i] * moleFractions[j] * Math.sqrt(componentParameters[i].E0 * componentParameters[j].E0);
            }
            parameter.B0 += moleFractions[i] * componentParameters[i].B0;
            parameter.a += moleFractions[i] * Math.cbrt(componentParameters[i].a);
            parameter.b += moleFractions[i] * Math.cbrt(componentParameters[i].b);
            parameter.c += moleFractions[i] * Math.cbrt(componentParameters[i].c);
            parameter.d += moleFractions[i] * Math.cbrt(componentParameters[i].d);
            parameter.alpha += moleFractions[i] * Math.cbrt(componentParameters[i].alpha);
            parameter.gamma += moleFractions[i] * Math.sqrt(componentParameters[i].gamma);
            parameter.Mg += moleFractions[i] * componentParameters[i].Mg;
        }
        parameter.B0 = Math.pow(parameter.B0, 3);
        parameter.a = Math.pow(parameter.a, 3);
        parameter.b = Math.pow(parameter.b, 3);
        parameter.c = Math.pow(parameter.c, 3);
        parameter.d = Math.pow(parameter.d, 3);
        parameter.alpha = Math.pow(parameter.alpha, 3);
        parameter.gamma = Math.pow(parameter.gamma, 2);
        return  parameter;
    }

    public static void main(String[] args) {
        double Tc = ComponentData.crititalTemperature_K[0];
        double denc = ComponentData.critialDensity_kmol_m3[0];
        double omega = ComponentData.acentricFactor[0];
        double Mg = ComponentData.molecularWeight_g_mol[0];
        double p = 101.325, t = 273.15;
        BWRS bwrs = new BWRS(BWRS.getParameter(Tc, denc, omega, Mg));
        double den = bwrs.getDensity(p, t) * ComponentData.molecularWeight_g_mol[0];
        System.out.println(den);
    }
}
