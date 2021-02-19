package cnic.pipesimulation;

import simiulatonutil.*;

import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.13 013 9:07:56
 */
public class SinglePipeSolver {
    public SinglePipeDiscretizer discretizer;
    public PropertyUpdater propertyUpdater;
    public SolverObserver observer = new DefaultObserver();
    public double deltaT;
    public double time = 0;
    public void simulate(double totalT) {
        // 温度方程相关参数
        double[] at, bt, ct, T, st;
        int nump = discretizer.nump;
        at = discretizer.at;
        bt = discretizer.bt;
        ct = discretizer.ct;
        T = discretizer.T;
        st = discretizer.st;
        // 动量方程相关参数
        double[] am, bm, cm, u, sm;
        int numv = discretizer.numv;
        am = discretizer.am;
        bm = discretizer.bm;
        cm = discretizer.cm;
        u = discretizer.u;
        sm = discretizer.sm;
        // 压力方程相关参数
        double[] ap, bp, cp, dp, sp;
        ap = discretizer.ap;
        bp = discretizer.bp;
        cp = discretizer.cp;
        dp = discretizer.dp;
        sp = discretizer.sp;

        observer.simulationStart(this);

        while (time < totalT) {
            observer.newTimestep(time);
            int itr = 0;
            // 内循环，求解当前时刻的速度压力
            while (discretizer.residual > 1e-5) {
                observer.innerLoopStart();
                // 计算物性
                propertyUpdater.propertyUpdate();
                observer.beforeDiscretizeMomentum();
                // 动量方程离散
                discretizer.momentumDiscretize();
                observer.afterDiscretizeMomentum();
                // 动量方程求解
                FunctionSolver.chase(am, bm, cm, u, 0, numv, sm);
                observer.beforeDiscretizePressure();
                // 压力方程离散
                discretizer.pressureDiscretize();
                observer.afterDiscretizePressure();
                if (Double.isNaN(discretizer.residual)) {
                    throw new RuntimeException("residual is NaN");
                }
                // 压力方程求解
                FunctionSolver.chase(ap, bp, cp, dp, 0, nump, sp);
                observer.beforeCorrect();
                // 修正速度压力
                discretizer.correct();
                observer.beforeDiscretizeTemperature();
                // 温度方程离散
                discretizer.temperatureDiscretize();
                observer.afterDiscretizeTemperature();
                // 温度方程求解
                FunctionSolver.chase(at, bt, ct, T, 0, nump, st);
                itr++;
                if (itr > 10) {
                    break;
                }
            }
            // 进入下一个时步
            discretizer.timeForward();
            time += deltaT;
            System.out.println("time = " + time);
        }
        observer.simulationEnd();
    }
    public static class Builder {
        SinglePipeDiscretizer discretizer;
        BWRS bwrs;
        CapacityCalculator capacity;
        Colebrook colebrook;
        PropertyUpdater propertyUpdater;
        BoundaryCondition boundary;
        double deltaX, deltaT;
        int numOfCells;
        double D, roughness, Area;
        double singleT, heatTransfer;
        double[] surroundT, heatTransfers;

        public Builder fluidComponents(String[] names, double[] moleFractions) {
            int len = names.length;
            int[] indices = new int[len];
            int n = ComponentData.formula.length;
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < n; j++) {
                    if (names[i] == ComponentData.formula[j]) {
                        indices[i] = j;
                        break;
                    }
                }
            }
            // 归一化
            double sum = 0;
            for (int i = 0; i < len; i++) {
                sum += moleFractions[i];
            }
            for (int i = 0; i < len; i++) {
                moleFractions[i] /= sum;
            }
            BWRS.ComponentParameter[] parameters = new BWRS.ComponentParameter[len];
            double[] crititalTemperature_K = ComponentData.crititalTemperature_K;
            double[] critialDensity_kmol_m3 = ComponentData.critialDensity_kmol_m3;
            double[] acentricFactor = ComponentData.acentricFactor;
            double[] molecularWeight_g_mol = ComponentData.molecularWeight_g_mol;
            for (int i = 0; i < len; i++) {
                int j = indices[i];
                parameters[i] = BWRS.getParameter(crititalTemperature_K[j], critialDensity_kmol_m3[j], acentricFactor[j], molecularWeight_g_mol[j]);
            }
            BWRS.ComponentParameter mixparameter = BWRS.getMixtureParameter(parameters, moleFractions);
            bwrs = new BWRS(mixparameter);

            double[] Mi = new double[len];
            for (int i = 0; i < len; i++) {
                Mi[i] = ComponentData.molecularWeight_g_mol[indices[i]];
            }
            double M = mixparameter.Mg;
            n = CapacityData.formula.length;
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < n; j++) {
                    if (names[i] == CapacityData.formula[j]) {
                        indices[i] = j;
                    }
                }
            }
            double B = 0, C = 0, D = 0, E = 0, F = 0, R = bwrs.R / M;
            for (int i = 0; i < len; i++) {
                B += moleFractions[i] * Mi[i] * CapacityData.B[indices[i]] / M;
                C += moleFractions[i] * Mi[i] * CapacityData.C[indices[i]] / M;
                D += moleFractions[i] * Mi[i] * CapacityData.D[indices[i]] / M;
                E += moleFractions[i] * Mi[i] * CapacityData.E[indices[i]] / M;
                F += moleFractions[i] * Mi[i] * CapacityData.F[indices[i]] / M;
            }
            double[] param = new double[]{B, C, D, E, F, R};
            capacity = new CapacityCalculator(param, bwrs);
            return this;
        }

        public Builder pipeSurrounding(double singleTemperature, double[] surroundingTemperature, double heatTransfer, double[] heatTransfers) {
            this.singleT = singleTemperature;
            this.surroundT = surroundingTemperature;
            this.heatTransfer = heatTransfer;
            this.heatTransfers = heatTransfers;
            return this;
        }

        public Builder discretize(double deltaT, double deltaX, int numOfCells) {
            this.deltaT = deltaT;
            this.deltaX = deltaX;
            this.numOfCells = numOfCells;
            return this;
        }

        public Builder pipeConfig(double D, double roughness) {
            this.D = D;
            this.roughness = roughness;
            colebrook = new Colebrook(roughness, D);
            Area = Math.PI * D * D / 4;
            return this;
        }

        public Builder boundary(BoundaryCondition boundary) {
            this.boundary = boundary;
            return this;
        }

        public SinglePipeSolver build() {
            SinglePipeSolver solver = new SinglePipeSolver();
            // 构建离散器
            discretizer = new SinglePipeDiscretizer();
            discretizer.initialArrays(numOfCells);
            discretizer.D = D;
            discretizer.Area = Area;
            discretizer.deltaX = deltaX;
            discretizer.deltaT = deltaT;
            discretizer.isFlowOut = boundary.isFlowOut;
            discretizer.isFlowIn = boundary.isFlowIn;
            discretizer.isPressIn = boundary.isPressIn;
            discretizer.isPressOut = boundary.isPressOut;
            discretizer.isTempIn = boundary.isTempIn;
            discretizer.isTempOut = boundary.isTempOut;
            Arrays.fill(discretizer.p, boundary.pressure);
            Arrays.fill(discretizer.p0, boundary.pressure);
            Arrays.fill(discretizer.T, boundary.temperature);
            Arrays.fill(discretizer.T0, boundary.temperature);
            if (surroundT == null) {
                Arrays.fill(discretizer.Ten, singleT);
            } else {
                System.arraycopy(surroundT, 0, discretizer.Ten, 0, numOfCells + 2);
            }
            if (heatTransfers == null) {
                Arrays.fill(discretizer.K, heatTransfer);
            } else {
                System.arraycopy(heatTransfers, 0, discretizer.K, 0, numOfCells + 2);
            }
            if (discretizer.isFlowIn) {
                discretizer.mass[1] = boundary.mass;
            }
            if (discretizer.isFlowOut) {
                discretizer.mass[numOfCells] = boundary.mass;
            }
            // 构建物性更新器
            propertyUpdater = new PropertyUpdater();
            propertyUpdater.discretizer = discretizer;
            propertyUpdater.bwrs = bwrs;
            propertyUpdater.capacity = capacity;
            propertyUpdater.colebrook = colebrook;
            // 构建求解器
            solver.discretizer = discretizer;
            solver.propertyUpdater = propertyUpdater;
            solver.deltaT = deltaT;
            // 初始条件
            Arrays.fill(discretizer.T, boundary.temperature);
            double den = bwrs.getDensity(boundary.pressure / 1000, boundary.temperature) * bwrs.parameter.Mg;
            Arrays.fill(discretizer.rho, den);
            Arrays.fill(discretizer.rho0, den);
            propertyUpdater.propertyUpdate();
            if (discretizer.isFlowIn) {
                discretizer.energy[1] = boundary.mass / Area / deltaX * discretizer.Cp[1] * boundary.temperature;
            }
            if (discretizer.isFlowOut) {
                discretizer.energy[numOfCells] = boundary.mass / Area / deltaX * discretizer.Cp[numOfCells] * boundary.temperature;
            }
            return solver;
        }
    }
}
