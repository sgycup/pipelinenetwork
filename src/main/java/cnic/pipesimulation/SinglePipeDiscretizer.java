package cnic.pipesimulation;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.08 008 16:01:46
 */
public class SinglePipeDiscretizer {
    // num为网格数，采用交错网格法与附加节点法离散方程，因此标量个数为nump=num+2，矢量个数为numv=num+1
    public int num, nump, numv;
    // 时间步长、空间步长、管径、横截面积
    public double deltaT, deltaX, D, Area;
    // 原始变量(压力、速度、密度、声速的平方、摩阻系数：按矢量处理)
    public double[] p, u, rho, c2, lambda;
    // 原始变量温度、定压比热、总换热系数、环境温度
    public double[] T, Cp, K, Ten;
    // 0表示上一时层的变量，dp表示压力修正量
    public double[] p0, u0, rho0, dp, T0;
    // 动量方程系统矩阵与源项
    public double[] am, bm, cm, sm;
    // 动量源(包括压缩机、阀等动力阻力元件，这里指这些元件产生的动力或阻力压差)
    public double[] deltaP;
    // 压力方程系数矩阵与源项
    public double[] ap, bp, cp, sp;
    // 压力修正与速度修正关系系数
    public double[] correctCoeff;
    // 质量源
    public double[] mass;// kg/s
    // 温度方程系数矩阵与源项
    public double[] at, bt, ct, st;
    // 能量源
    public double[] energy;// W/m3
    // 出入口定流标志
    public boolean isFlowOut, isFlowIn;
    // 出入口定压标志
    public boolean isPressOut, isPressIn;
    // 出入口定温标志
    public boolean isTempOut, isTempIn;
    // 残差
    public double residual = Double.MAX_VALUE;

    public void initialArrays(int num) {
        this.num = num;
        this.numv = num + 1;
        this.nump = num + 2;
        p = new double[nump];
        u = new double[numv];
        rho = new double[nump];
        c2 = new double[nump];
        lambda = new double[numv];
        T = new double[nump];
        Ten = new double[nump];
        energy = new double[nump];
        Cp = new double[nump];
        K = new double[nump];
        p0 = new double[nump];
        u0 = new double[numv];
        rho0 = new double[nump];
        dp = new double[nump];
        T0 = new double[nump];
        am = new double[numv - 1];
        bm = new double[numv];
        cm = new double[numv - 1];
        sm = new double[numv];
        ap = new double[nump - 1];
        bp = new double[nump];
        cp = new double[nump - 1];
        sp = new double[nump];
        correctCoeff = new double[numv];
        deltaP = new double[numv];
        mass = new double[nump];
        at = new double[nump - 1];
        bt = new double[nump];
        ct = new double[nump - 1];
        st = new double[nump];
    }

    public void momentumDiscretize() {
        // 中间节点
        for (int i = 1, n = numv - 1; i < n; i++) {
            // 非稳态项
            double rhoB = (rho[i] + rho[i + 1]) / 2;
            bm[i] = rhoB / deltaT;
            // 对流项
            double ub = u[i] / deltaX;
            if (ub >= 0) {
                bm[i] += rho[i + 1] * ub;
                am[i - 1] = -rho[i] * ub;
                cm[i] = 0;
            } else {
                bm[i] -= rho[i] * ub;
                cm[i] = rho[i + 1] * ub;
                am[i - 1] = 0;
            }
            // 源项
            // 非稳态项
            double rho0B = (rho0[i] + rho0[i + 1]) / 2;
            sm[i] = rho0B * u0[i] / deltaT;
            // 压力项
            sm[i] += -(p[i + 1] - p[i] - deltaP[i]) / deltaX;
            // 摩阻
            bm[i] += lambda[i] * rhoB * Math.abs(u[i]) / 2 / D;
            // 质量源
            if (mass[i] * ub >= 0) {
                sm[i] += mass[i] / Area / deltaX;
            }
            if (mass[i + 1] * ub <= 0) {
                sm[i] += -mass[i + 1] / Area / deltaX;
            }
            // 修正系数
            correctCoeff[i] = -1 / deltaX / bm[i];
        }
        // 入口节点
        // 非稳态项
        double rhoB = (rho[0] + rho[1]) / 2;
        bm[0] = rhoB / deltaT;
        // 对流项
        double ub = u[0] / deltaX;
        bm[0] -= rho[0] * ub;
        cm[0] = rho[1] * ub;
        // 源项
        // 非稳态项
        double rho0B = (rho0[0] + rho0[1]) / 2;
        sm[0] = rho0B * u0[0] / deltaT;
        // 压力项
        sm[0] += -(p[1] - p[0]) / deltaX;
        // 摩阻
        bm[0] += lambda[0] * rhoB * Math.abs(u[0]) / 2 / D;
        // 修正系数
        correctCoeff[0] = -1 / deltaX / bm[0];
        // 质量源
        if (isFlowIn) {
            // 将第一个节点置零
            bm[0] = 1;
            cm[0] = 0;
            sm[0] = 0;
            correctCoeff[0] = 0;
        }
        // 出口节点
        // 非稳态项
        int out = numv - 1;
        rhoB = (rho[out] + rho[out + 1]) / 2;
        bm[out] = rhoB / deltaT;
        // 对流项
        ub = u[out] / deltaX;
        bm[out] += rho[out + 1] * ub;
        am[out - 1] = -rho[out] * ub;
        // 源项
        // 非稳态项
        rho0B = (rho0[out] + rho0[out + 1]) / 2;
        sm[out] = rho0B * u0[out] / deltaT;
        // 压力项
        sm[out] += -(p[out + 1] - p[out]) / deltaX;
        // 摩阻
        bm[out] += lambda[out] * rhoB * Math.abs(u[out]) / 2 / D;
        // 修正系数
        correctCoeff[out] = -1 / deltaX / bm[out];
        // 质量源
        if (isFlowOut) {
            // 将最后一个节点置零
            bm[out] = 1;
            am[out - 1] = 0;
            sm[out] = 0;
            correctCoeff[out] = 0;
        }
    }

    public void pressureDiscretize() {
        bp[0] = deltaX / c2[0] / deltaT;
        for (int i = 0, ix; i < numv; i++) {
            ix = i + 1;
            double rhoB = u[i] >= 0 ? rho[i] : rho[ix];
            double tmp = rhoB * correctCoeff[i];
            ap[i] = tmp;
            cp[i] = tmp;
            bp[i] -= tmp;
            bp[ix] = -tmp;
            bp[ix] += deltaX / c2[ix] / deltaT;
        }

        for (int i = 1; i < numv; i++) {
            double rhoBL = u[i - 1] >= 0 ? rho[i - 1] : rho[i];
            double rhoBR = u[i] >= 0 ? rho[i] : rho[i + 1];
            sp[i] = rhoBL * u[i - 1] - rhoBR * u[i];
            sp[i] += deltaX * (rho0[i] - rho[i]) / deltaT;
            if (mass[i] != 0) {
                sp[i] += mass[i] / Area;
            }
        }

        if (isPressIn) {
            bp[0] = 1;
            cp[0] = 0;
            sp[0] = 0;
        }
        if (isPressOut) {
            bp[numv] = 1;
            ap[num] = 0;
            sp[numv] = 0;
        }
        if (isFlowIn) {
            bp[0] = 1;
            cp[0] = -1;
            sp[0] = 0;
        }
        if (isFlowOut) {
            bp[numv] = 1;
            ap[num] = -1;
            sp[numv] = 0;
        }

        residual = 0;
        for (int i = 1; i < numv; i++) {
            residual += sp[i] * sp[i];
        }
        residual = Math.sqrt(residual);
    }

    public void correct() {
        for (int i = 0; i < nump; i++) {
            p[i] += dp[i];
        }
        for (int i = 0; i < nump; i++) {
            rho[i] += dp[i] / c2[i];
            if (rho[i] < 0) {
                throw new RuntimeException("rho < 0 (rho = " + rho[i] + ", i = " + i + ")");
            }
        }
        for (int i = 0; i < numv; i++) {
            u[i] += (dp[i + 1] - dp[i]) * correctCoeff[i];
        }
    }

    public void temperatureDiscretize() {
        // 质量源引起的热源变化
        if (isFlowIn) {
            energy[1] = mass[1] / Area / deltaX * Cp[1] * T[1];
        }
        if (isFlowOut) {
            energy[num] = mass[num] / Area / deltaX * Cp[num] * T[num];
        }
        // 温度方程离散
        for (int i = 1; i < numv; i++) {
            // 非稳态项
            bt[i] = rho[i] * Cp[i] / deltaT;
            // 网格左侧
            if (u[i - 1] >= 0) {
                at[i - 1] = -rho[i - 1] * Cp[i] * u[i - 1] / deltaX;
            } else {
                bt[i] -= rho[i] * Cp[i] * u[i - 1] / deltaX;
                at[i - 1] = 0;
            }
            // 网格右侧
            if (u[i] >= 0) {
                bt[i] += rho[i] * Cp[i] * u[i] / deltaX;
                ct[i] = 0;
            } else {
                ct[i] = rho[i + 1] * Cp[i] * u[i] / deltaX;
            }

            // 源项
            // 非稳态项
            st[i] = rho[i] * Cp[i] * T0[i] / deltaT;
            // 与环境换热，采用隐式离散
            double KPiD = K[i] * Math.PI * D / Area;
            bt[i] += KPiD;
            st[i] += KPiD * Ten[i];
            // 推动功
            double vc = (u[i - 1] + u[i]) / 2;
            st[i] += vc * (p[i + 1] - p[i - 1]) / deltaX / 2;
            // 摩阻做功
            st[i] += lambda[i] * rho[i] * Math.abs(Math.pow(vc, 3)) / 2;
            // 热源
            st[i] += energy[i];
        }
        // 边界条件
        if (isTempIn) {
            bt[0] = 1;
            ct[0] = 0;
            st[0] = T[0];
        } else {
            bt[0] = 1;
            ct[0] = -1;
            st[0] = 0;
        }
        if (isTempOut) {
            bt[numv] = 1;
            at[num] = 0;
            st[numv] = T[numv];
        } else {
            bt[numv] = 1;
            at[num] = -1;
            st[numv] = 0;
        }
    }

    public void timeForward() {
        System.arraycopy(u, 0, u0, 0, numv);
        System.arraycopy(p, 0, p0, 0, nump);
        System.arraycopy(rho, 0, rho0, 0, nump);
        System.arraycopy(T, 0, T0, 0, nump);
        residual = Double.MAX_VALUE;
    }
}
