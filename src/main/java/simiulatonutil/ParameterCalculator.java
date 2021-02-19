package simiulatonutil;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.12 012 9:57:47
 */
public class ParameterCalculator {
    private ParameterCalculator(){}

    /**
     * @param density 密度
     * @param u 速度
     * @param viscosity 动力粘度
     * @param len 参考长度
     * @return 雷诺数
     */
    public static double Reynolds(double density, double u, double viscosity, double len) {
        return density * u * len / viscosity;
    }

    /**
     * @param u 速度
     * @param viscosity 运动粘度
     * @param len 参考长度
     * @return 雷诺数
     */
    public static double Reynolds(double u, double viscosity, double len) {
        return u * len / viscosity;
    }

    /**
     * @param density 密度，kg/m3
     * @param temperature 温度，K
     * @return 粘度，Pa*s
     */
    public static double viscosity(double density, double temperature) {
        double delta = density / 1.293;
        double x = 2.57 + 0.2781 * delta + 1063.6 / temperature;
        double y = 1.11 + 0.04 * x;
        double C = 2.415 * (7.77 + 0.1844 * delta) * Math.pow(temperature, 1.5);
        C /= 122.4 + 377.58 * delta + 1.8 * temperature;
        C *= 1e-4;
        return C * Math.exp(x * Math.pow(density / 1000, y)) / 1000;
    }

    public static double lambdaColebrook(Colebrook f, double Re) {
        if (Re < 1e-5) {
            return 0;
        } else if (Re < 2000) {
            return 64 / Re;
        } else if (Re < 4000) {
            return 0.0025 * Math.cbrt(Re);
        } else {
            f.setRe(Re);
            return FunctionSolver.secantSolver(f, 0, 0.1, 1e-6);
        }
    }

    /**
     * @Author Shi Guoyun
     * @Description 将质量分数转化为摩尔分数
     * @Date 21:33:21 2021.01.25 025
     *
     * @param names 组分名称
     * @param fractions 组分质量分数
     * @return 摩尔分数
     */
    public static double[] weight2mole(String[] names, double[] fractions) {
        int len = names.length;
        double[] results = new double[len];
        double[] Mg = new double[len];
        String[] formula = ComponentData.formula;
        int n = formula.length;
        for (int i = 0; i < len; i++) {
            for (int j = 0; j < n; j++) {
                if (names[i] == formula[j]) {
                    Mg[i] = ComponentData.molecularWeight_g_mol[j];
                    break;
                }
            }
        }
        double sum = 0;
        for (int i = 0; i < len; i++) {
            sum += fractions[i] / Mg[i];
            results[i] = fractions[i] / Mg[i];
        }
        for (int i = 0; i < len; i++) {
            results[i] /= sum;
        }
        return results;
    }

    public static double[] normalize(double[] fractions) {
        int len = fractions.length;
        double[] results = new double[len];
        double sum = 0;
        for (int i = 0; i < len; i++) {
            sum += fractions[i];
        }
        for (int i = 0; i < len; i++) {
            results[i] = fractions[i] / sum;
        }
        return results;
    }
}
