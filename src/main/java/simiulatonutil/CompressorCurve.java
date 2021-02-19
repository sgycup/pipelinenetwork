package simiulatonutil;

/**
 * @Description: 压缩机不同转速下特性曲线(压比 - 流量)
 * epsilon^2 = A - B * Q^2
 * Q的单位为：am3/h，工况下立方米每小时。epsilon = P2 / P1，出口压力比入口压力
 *
 * 还给出了这几条特性曲线回归时的最大最小流量。
 *
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.27 027 15:58:26
 */
public class CompressorCurve {
    public static final double[] RPM = new double[]{4412, 4751, 5430, 6109, 6788, 7127};
    public static final double[] A = new double[]{2.18863751334086, 2.0425540828574245, 1.7975747731299418, 1.6237249123232145, 1.4501641634749491, 1.374954131849213};
    public static final double[] B = new double[]{1.3733682435635247E-8, 1.0005405229529623E-8, 5.38635444032285E-9, 3.4461421073442903E-9, 1.8255814107563319E-9, 1.3072767467791275E-9};
    public static final double[] FLOWRATE_MIN = new double[]{3741, 4044, 4648, 5369, 6235, 6719};
    public static final double[] FLOWRATE_MAX = new double[]{7441, 8065, 9314, 10619, 12000, 12685};
}
