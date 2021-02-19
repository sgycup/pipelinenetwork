package simiulatonutil;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.12 012 16:01:47
 */
public class ComponentData {
    public final static String[] formula = {"CH4", "C2H4", "C2H6", "C3H6", "C3H8", "i-C4H10", "n-C4H10", "i-C5H12",
            "n-C5H12", "C6H14", "C7H16", "C8H18", "C9H20", "C10H22", "C11H24", "N2", "CO2", "H2S", "H2", "H20"};
    public final static double[] crititalPressure_MPa = {4.604, 5.03, 4.88, 4.63, 4.25, 3.648, 3.797, 3.374,
            3.369, 3.012, 2.736, 2.487, 2.31, 2.09, 1.95, 3.394, 7.376, 9, 13.2, 22.048};
    public final static double[] crititalTemperature_K = {190.69, 283.05, 305.38, 365.04, 369.89, 408.13, 425.18, 460.37,
            469.49, 507.28, 540.28, 568.58, 594.57, 617.54, 639.99, 126.15, 304.09, 373.39, 47.05, 647.3};
    public final static double[] critialDensity_kmol_m3 = {10.05, 8.0653, 6.7566, 5.5248, 4.9994, 3.8012, 3.9213, 3.2469,
            3.2149, 2.7167, 2.3467, 2.0568, 1.8421, 1.6611, 1.5154, 11.099, 10.638, 10.526, 20, 17.857};
    public final static double[] acentricFactor = {0.013, 0.101, 0.1018, 0.15, 0.157, 0.183, 0.197, 0.226, 0.252,
            0.302, 0.353, 0.412, 0.475, 0.54, 0.6, 0.035, 0.21, 0.105, 0, 0.344};
    public final static double[] molecularWeight_g_mol = {16.042, 28.05, 30.068, 42.08, 44.094, 58.12, 58.12, 72.146,
            72.146, 86.172, 100.198, 114.224, 128.25, 142.276, 156.3, 28.016, 44.01, 34.076, 2, 18.015};
}
