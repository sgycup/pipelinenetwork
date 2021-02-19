package simiulatonutil;

import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.08 008 10:18:13
 */
public class FunctionSolver {
    /**
     * 追赶法求解三对角方程组
     *
     * @param a     下副对角线
     * @param b     主对角线
     * @param c     上副对角线
     * @param x     未知变量
     * @param start 子矩阵对角线起始索引
     * @param end   子矩阵对角线终止索引
     * @param d     源项
     * @param ed    额外源项，方法结束时变为对应的解
     * @return 源项对应的解，当对角线第二个元素为0时直接返回初始值。采用第二元素判断的原因是第一个元素有可能为边界条件，使得其不会为0，不能反应实际情况
     */
    public static double[] chase(double[] a, double[] b, double[] c, double[] x, int start, int end, double[] d, double[]... ed) {
        if (b[start + 1] == 0) {
            return x;
        }
        // 分解
        for (int i = start + 1; i < end; i++) {
            a[i - 1] /= b[i - 1];
            b[i] -= c[i - 1] * a[i - 1];
        }
        // 解下三角
        x[start] = d[start];
        for (int i = start + 1; i < end; i++) {
            x[i] = d[i] - a[i - 1] * x[i - 1];
        }
        for (double[] dj : ed) {
            for (int i = start + 1; i < end; i++) {
                dj[i] = dj[i] - a[i - 1] * dj[i - 1];
            }
        }
        // 解上三角
        x[end - 1] /= b[end - 1];
        for (int i = end - 2; i >= start; i--) {
            x[i] = (x[i] - c[i] * x[i + 1]) / b[i];
        }
        for (double[] dj : ed) {
            dj[end - 1] /= b[end - 1];
            for (int i = end - 2; i >= start; i--) {
                dj[i] = (dj[i] - c[i] * dj[i + 1]) / b[i];
            }
        }
        return x;
    }

    public static double bisectionSolver(UnivariateFunction f, double min, double max, double absoluteAccuracy) {
        double m;
        double fm;
        double fmin;
        while (true) {
            m = (min + max) / 2;
            fmin = f.value(min);
            fm = f.value(m);

            if (fm * fmin > 0) {
                min = m;
            } else {
                max = m;
            }

            if (Math.abs(max - min) <= absoluteAccuracy) {
                m = (min + max) / 2;
                return m;
            }
        }
    }

    public static double secantSolver(UnivariateFunction f, double min, double max, double absoluteAccuracy) {
        double x0 = min;
        double x1 = max;
        double f0 = f.value(x0);
        double f1 = f.value(x1);

        if (f0 == 0.0) {
            return x0;
        }
        if (f1 == 0.0) {
            return x1;
        }


        while (true) {

            final double x = x1 - ((f1 * (x1 - x0)) / (f1 - f0));
            final double fx = f.value(x);

            if (fx == 0.0) {
                return x;
            }

            x0 = x1;
            f0 = f1;
            x1 = x;
            f1 = fx;

            if (Math.abs(x1 - x0) < absoluteAccuracy) {
                return x1;
            }
        }
    }

    public static void sum(double a, double[] x, double b, double[] y, double[] sum) {
        int n = x.length;
        for (int i = 0; i < n; i++) {
            sum[i] = a * x[i] + b * y[i];
        }
    }

    public static double dot(double[] x, double[] y) {
        double prod = 0;
        int n = x.length;
        for (int i = 0; i < n; i++) {
            prod += x[i] * y[i];
        }
        return prod;
    }

    public static double norm(double[] x) {
        double sum = 0;
        for (double v : x) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }

    public static class BiCGStab {
        double[] r0, r, v, p, s, t;
        double rho0, rho, alpha, beta, omega;
        double threshold;
        int maxitr;

        public BiCGStab(int len, int maxitr, double threshold) {
            this.maxitr = maxitr;
            this.threshold = threshold;
            r0 = new double[len];
            r = new double[len];
            v = new double[len];
            p = new double[len];
            s = new double[len];
            t = new double[len];
        }

        public void solve(Matrix A, double[] x, double[] b) {
            A.preprocess(b);
            A.multipy(x, r0);
            sum(1, b, -1, r0, r0);
            if (norm(r0) < threshold) {
                return;
            }
            int n = r0.length;
            System.arraycopy(r0, 0, r, 0, n);
            rho0 = alpha = omega = 1;
            Arrays.fill(v, 0);
            Arrays.fill(p, 0);
            for (int i = 0; i < maxitr; i++) {
                rho = dot(r, r0);
                beta = rho / rho0 * alpha / omega;
                sum(1, p, -omega, v, p);
                sum(1, r, beta, p, p);
                A.multipy(p, v);
                alpha = rho / dot(r0, v);
                sum(1, r, -alpha, v, s);
                if (norm(s) < threshold) {
                    sum(1, x, alpha, p, x);
                    break;
                }
                A.multipy(s, t);
                omega = dot(t, s) / dot(s, s);
                sum(1, x, alpha, p, x);
                sum(1, x, omega, s, x);
                sum(1, s, -omega, t, r);
                rho0 = rho;
            }
        }
    }

    public static void main(String[] args) {
        int n = 3;
        double[] a = new double[n - 1];
        double[] b = new double[n];
        double[] c = new double[n - 1];
        double[] d = new double[n];
        double[] ed = new double[n];
        double[] x = new double[n];
        Arrays.fill(a, -1);
        Arrays.fill(c, -1);
        Arrays.fill(b, 4);
        b[0] = b[n - 1] = 3;
        Arrays.fill(d, 2);
        Arrays.fill(ed, 2);
        FunctionSolver.chase(a, b, c, x, 0, n, d, ed);
        System.out.println(Arrays.toString(x));
        System.out.println(Arrays.toString(ed));
    }
}
