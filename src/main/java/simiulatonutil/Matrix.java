package simiulatonutil;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.08 008 10:17:11
 */

public interface Matrix {

    /**
     * 预处理
     *
     * @param s 源项
     */
    default void preprocess(double[] s){}

    /**
     * y = A * x
     *
     * @param x 向量
     * @param y 向量
     */
    void multipy(double[] x, double[] y);
}
