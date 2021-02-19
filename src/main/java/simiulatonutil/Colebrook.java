package simiulatonutil;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.12 012 19:47:24
 */
public class Colebrook implements UnivariateFunction {

    // 绝对粗糙度与管径
    private final double Ke, D;

    // 雷诺数
    private double Re;

    public Colebrook(double ke, double d) {
        Ke = ke;
        D = d;
    }

    @Override
    public double value(double x) {
        double tmp = 1 / Math.sqrt(x);
        return tmp + 2 * Math.log10(Ke / 3.7 / D + 2.51 / Re * tmp);
    }

    public void setRe(double re) {
        Re = re;
    }

}
