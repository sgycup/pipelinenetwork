package cnic.pipesimulation;

import simiulatonutil.BWRS;
import simiulatonutil.CSRMatrix;
import simiulatonutil.FunctionSolver;
import simiulatonutil.PropertyUpdater;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.10 010 14:59:52
 */
public class NetworkSolver {
    public double[] rho, u, dp, T, volume;// 边界点体积为0
    public double[] rho0, T0;
    public Node[] nodes;
    public Edge[] edges;
    public double residual = Double.MAX_VALUE;
    public double deltaT;
    public double time = 0;
    public double EnergyStart = 0;

    public void simulate(double timeTotal) {
        int numEdge = edges.length;
        int numNode = nodes.length;

        CSRMatrix csrv = createVelocityMatrix();
        double[] sourcev = new double[numEdge * 2];

        CSRMatrix csrp = createPressureMatrix();
        double[] sourcep = new double[numNode];

        CSRMatrix csrt = createPressureMatrix();
        double[] sourcet = new double[numNode];

        FunctionSolver.BiCGStab biCGStabV = new FunctionSolver.BiCGStab(numEdge * 2, 1000, 1e-6);
        FunctionSolver.BiCGStab biCGStabP = new FunctionSolver.BiCGStab(numNode, 1000, 1e-6);
        FunctionSolver.BiCGStab biCGStabT = new FunctionSolver.BiCGStab(numNode, 1000, 1e-6);

        for (int i = 0; i < numNode; i++) {
            Node node = nodes[i];
            boolean isSource = node.isSource.get(0);
            Edge edge = edges[node.edgeID.get(0)];
            double[] den = edge.pipeDiscretizer.rho;
            double[] temp = edge.pipeDiscretizer.T;
            int num = edge.pipeDiscretizer.num;
            if (isSource) {
                rho[i] = rho0[i] = den[0];
                T[i] = T0[i] = temp[0];
            } else {
                rho[i] = rho0[i] = den[num + 1];
                T[i] = T0[i] = temp[num + 1];
            }
        }

        while (time < timeTotal) {
            int itr = 0;
            // 内循环，求解当前时刻的速度压力
            while (residual > 1e-5) {
                // 计算物性
                for (int i = 0; i < numEdge; i++) {
                    Edge edge = edges[i];
                    edge.propertyUpdater.propertyUpdate();
                }
                // 离散各管道内部节点动量方程
                for (int i = 0; i < numEdge; i++) {
                    double[][] freev = edges[i].freev;
                    Arrays.fill(freev[0], 0);
                    Arrays.fill(freev[1], 0);
                }
                for (int i = 0; i < numEdge; i++) {
                    Edge edge = edges[i];
                    SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                    int num = singlePipe.num;
                    singlePipe.momentumDiscretize();
                    double[][] freev = edge.freev;
                    freev[0][1] = -singlePipe.am[0];
                    freev[1][num - 1] = -singlePipe.cm[num - 1];
                    FunctionSolver.chase(singlePipe.am, singlePipe.bm, singlePipe.cm, singlePipe.u, 1, num, singlePipe.sm, freev);
                }
                // 离散网络节点动量方程
                csrv.clearData();
                for (int i = 0; i < numNode; i++) {
                    Node node = nodes[i];
                    ArrayList<Integer> edgeID = node.edgeID;
                    int n = edgeID.size();
                    // 边界节点处理
                    if (n == 1) {
                        boundaryMomentum(csrv, sourcev, node);
                        continue;
                    }
                    // 每个节点相连的边
                    nodeMomentum(csrv, sourcev, i, node);
                }

                // 求解网络节点动量方程
                biCGStabV.solve(csrv, u, sourcev);
                // 联立得到整体网络速度
                for (int i = 0; i < numEdge; i++) {
                    Edge edge = edges[i];
                    SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                    int num = singlePipe.num;
                    double x0 = u[2 * i];
                    double xn = u[2 * i + 1];
                    double[] pipeu = singlePipe.u;
                    pipeu[0] = x0;
                    pipeu[num] = xn;
                    double[][] freev = edge.freev;
                    double[] f0 = freev[0];
                    double[] fn = freev[1];
                    for (int j = 1; j < num; j++) {
                        pipeu[j] += x0 * f0[j] + xn * fn[j];
                    }
                }

                residual = 0;
                // 求解各管道内部节点的压力
                for (int i = 0; i < numEdge; i++) {
                    double[][] freep = edges[i].freep;
                    Arrays.fill(freep[0], 0);
                    Arrays.fill(freep[1], 0);
                }
                for (int i = 0; i < numEdge; i++) {
                    Edge edge = edges[i];
                    SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                    int num = singlePipe.num;
                    singlePipe.pressureDiscretize();
                    residual += singlePipe.residual * singlePipe.residual;
                    double[][] freep = edge.freep;
                    freep[0][1] = -singlePipe.ap[0];
                    freep[1][num] = -singlePipe.cp[num];
                    FunctionSolver.chase(singlePipe.ap, singlePipe.bp, singlePipe.cp, singlePipe.dp, 1, num + 1, singlePipe.sp, freep);
                }

                // 离散网络节点压力方程
                csrp.clearData();
                for (int i = 0; i < numNode; i++) {
                    Node node = nodes[i];
                    // 节点定压
                    if (node.pressure != 0) {
                        csrp.setData(i, i, 1);
                        sourcep[i] = 0;
                        continue;
                    }
                    ArrayList<Integer> edgeID = node.edgeID;
                    int n = edgeID.size();
                    // 边界节点处理
                    if (n == 1) {
                        boundaryPressure(csrp, sourcep, i, node);
                        continue;
                    }
                    // 节点压力
                    nodePressure(csrp, sourcep, i, node);
                }
                residual = Math.sqrt(residual);

                // 求解网络节点压力方程
                biCGStabP.solve(csrp, dp, sourcep);
                // 联立得到整体网络压力
                for (int i = 0; i < numEdge; i++) {
                    Edge edge = edges[i];
                    SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                    int numv = singlePipe.numv;
                    double x0 = dp[edge.source];
                    double xn = dp[edge.target];
                    double[] pipep = singlePipe.dp;
                    pipep[0] = x0;
                    pipep[numv] = xn;
                    double[][] freep = edge.freep;
                    double[] f0 = freep[0];
                    double[] fn = freep[1];
                    for (int j = 1; j < numv; j++) {
                        pipep[j] += x0 * f0[j] + xn * fn[j];
                    }
                }

                // 更新变量
                for (int i = 0; i < numEdge; i++) {
                    Edge edge = edges[i];
                    SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                    singlePipe.correct();
                }
                for (int i = 0; i < numNode; i++) {
                    Node node = nodes[i];
                    boolean isSource = node.isSource.get(0);
                    Edge edge = edges[node.edgeID.get(0)];
                    double[] den = edge.pipeDiscretizer.rho;
                    int num = edge.pipeDiscretizer.num;
                    if (isSource) {
                        rho[i] = den[0];
                    } else {
                        rho[i] = den[num + 1];
                    }
                }

                if (time >= EnergyStart) {
                    // 求解各管道内部节点的温度
                    for (int i = 0; i < numEdge; i++) {
                        double[][] freet = edges[i].freet;
                        Arrays.fill(freet[0], 0);
                        Arrays.fill(freet[1], 0);
                    }
                    for (int i = 0; i < numEdge; i++) {
                        Edge edge = edges[i];
                        SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                        int num = singlePipe.num;
                        singlePipe.temperatureDiscretize();
                        double[][] freet = edge.freet;
                        freet[0][1] = -singlePipe.at[0];
                        freet[1][num] = -singlePipe.ct[num];
                        FunctionSolver.chase(singlePipe.at, singlePipe.bt, singlePipe.ct, singlePipe.T, 1, num + 1, singlePipe.st, freet);
                    }
                    // 离散网络节点温度方程
                    csrt.clearData();
                    for (int i = 0; i < numNode; i++) {
                        Node node = nodes[i];
                        // 节点定压
                        if (node.temperature != 0) {
                            csrt.setData(i, i, 1);
                            sourcet[i] = node.temperature;
                            continue;
                        }
                        ArrayList<Integer> edgeID = node.edgeID;
                        int n = edgeID.size();
                        // 边界点处理
                        if (n == 1) {
                            boundaryTemperature(csrt, sourcet, i, node);
                            continue;
                        }
                        // 节点温度
                        nodeTemperature(csrt, sourcet, i, node);
                    }
                    // 求解网络节点温度方程
                    biCGStabT.solve(csrt, T, sourcet);
                    // 联立得到整体网络温度
                    for (int i = 0; i < numEdge; i++) {
                        Edge edge = edges[i];
                        SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                        int numv = singlePipe.numv;
                        double x0 = T[edge.source];
                        double xn = T[edge.target];
                        double[] pipet = singlePipe.T;
                        pipet[0] = x0;
                        pipet[numv] = xn;
                        double[][] freet = edge.freet;
                        double[] f0 = freet[0];
                        double[] fn = freet[1];
                        for (int j = 1; j < numv; j++) {
                            pipet[j] += x0 * f0[j] + xn * fn[j];
                        }
                    }
                }

                itr++;
                if (itr > 10) {
                    break;
                }
            }
            // 更新时步
            for (int i = 0; i < numEdge; i++) {
                Edge edge = edges[i];
                SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
                singlePipe.timeForward();
            }
            System.arraycopy(rho, 0, rho0, 0, numNode);
            System.arraycopy(T, 0, T0, 0, numNode);
            time += deltaT;
            residual = Double.MAX_VALUE;
            System.out.println("time = " + time);
        }
    }

    private void nodeTemperature(CSRMatrix csrt, double[] sourcet, int i, Node node) {
        ArrayList<Boolean> isSource = node.isSource;
        ArrayList<Integer> edgeID = node.edgeID;
        int n = edgeID.size();
        double deltaT = 1, rhocp = 0, cp = 0, V = volume[i];
        sourcet[i] = 0;
        double aveK = 0, aveA = 0, aveD = 0, Ten = 0;
        // 对流项
        for (int j = 0; j < n; j++) {
            int id = edgeID.get(j);
            Edge edge = edges[id];
            SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
            deltaT = singlePipe.deltaT;
            int num = singlePipe.num;
            double area = singlePipe.Area;
            aveA += area;
            aveD += singlePipe.D;
            double[][] freet = edge.freet;
            if (isSource.get(j)) {
                aveK += singlePipe.K[0];
                Ten += singlePipe.Ten[0];
                cp = singlePipe.Cp[0];
                rhocp = singlePipe.rho[0] * singlePipe.Cp[0];
                double ub = singlePipe.u[0];
                if (ub >= 0) {
                    csrt.addData(i, i, rhocp * ub * area / V);
                } else {
                    double c = singlePipe.rho[1] * singlePipe.Cp[1] * ub * area / V;
                    double f0 = freet[0][1];
                    double fn = freet[1][1];
                    double x1 = singlePipe.T[1];
                    csrt.addData(i, i, c * f0);
                    csrt.addData(i, edge.target, c * fn);
                    sourcet[i] -= c * x1;
                }
            } else {
                aveK += singlePipe.K[num + 1];
                Ten += singlePipe.Ten[num + 1];
                cp = singlePipe.Cp[num + 1];
                rhocp = singlePipe.rho[num + 1] * singlePipe.Cp[num + 1];
                double ub = singlePipe.u[num];
                if (ub >= 0) {
                    double a = -singlePipe.rho[num] * singlePipe.Cp[num] * ub * area / V;
                    double f0 = freet[0][num];
                    double fn = freet[1][num];
                    double xn = singlePipe.T[num];
                    csrt.addData(i, i, a * fn);
                    csrt.addData(i, edge.source, a * f0);
                    sourcet[i] -= a * xn;
                } else {
                    csrt.addData(i, i, -rhocp * ub * area / V);
                }
            }
        }
        // 系数矩阵与源项中的非稳态项
        csrt.addData(i, i, rhocp / deltaT);
        sourcet[i] += rhocp * T0[i] / deltaT;
        // 与环境换热项
        aveD /= n;
        aveA = Math.PI * aveD * aveD / 4;
        aveK /= n;
        Ten /= n;
        double b = Math.PI * aveK * aveD / aveA;
        csrt.addData(i, i, b);
        sourcet[i] += b * Ten;
        // 热源
        sourcet[i] += node.mass / V * cp * T[i];
    }

    private void boundaryTemperature(CSRMatrix csrt, double[] sourcet, int i, Node node) {
        ArrayList<Boolean> isSource = node.isSource;
        ArrayList<Integer> edgeID = node.edgeID;
        int id = edgeID.get(0);
        Edge edge = edges[id];
        SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
        double[][] freet = edge.freet;
        int num = singlePipe.num;
        if (isSource.get(0)) {
            sourcet[i] = singlePipe.st[0];
            double b = singlePipe.bt[0];
            double c = singlePipe.ct[0];
            double f0 = freet[0][1];
            double fn = freet[1][1];
            double x1 = singlePipe.T[1];
            csrt.setData(i, i, b + c * f0);
            csrt.setData(i, edge.target, c * fn);
            sourcet[i] -= c * x1;
        } else {
            sourcet[i] = singlePipe.st[num + 1];
            double b = singlePipe.bt[num + 1];
            double a = singlePipe.at[num];
            double f0 = freet[0][num];
            double fn = freet[1][num];
            double xn = singlePipe.T[num];
            csrt.setData(i, i, b + a * fn);
            csrt.setData(i, edge.source, a * f0);
            sourcet[i] -= a * xn;
        }
    }

    private void nodePressure(CSRMatrix csrp, double[] sourcep, int i, Node node) {
        ArrayList<Boolean> isSource = node.isSource;
        ArrayList<Integer> edgeID = node.edgeID;
        int n = edgeID.size();
        double deltaT = 1, c2 = 1;
        sourcep[i] = 0;
        double residualtmp = 0;
        for (int j = 0; j < n; j++) {
            int id = edgeID.get(j);
            Edge edge = edges[id];
            SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
            deltaT = singlePipe.deltaT;
            int num = singlePipe.num;
            double area = singlePipe.Area;
            double[][] freep = edge.freep;
            if (isSource.get(j)) {
                c2 = singlePipe.c2[0];
                rho[i] = singlePipe.rho[0];
                rho0[i] = singlePipe.rho0[0];
                double tmp = singlePipe.ap[0] * area;
                csrp.addData(i, edge.target, tmp * freep[1][1]);
                csrp.addData(i, i, tmp * (freep[0][1] - 1));
                double rhoB = singlePipe.u[0] >= 0 ? singlePipe.rho[0] : singlePipe.rho[1];
                double flux = -rhoB * singlePipe.u[0] * area;
                sourcep[i] += flux;
                residualtmp += flux;
                double x1 = singlePipe.dp[1];
                sourcep[i] += -tmp * x1;
            } else {
                c2 = singlePipe.c2[num + 1];
                rho[i] = singlePipe.rho[num + 1];
                rho0[i] = singlePipe.rho0[num + 1];
                double tmp = -singlePipe.cp[num] * area;
                csrp.addData(i, edge.source, -tmp * freep[0][num]);
                csrp.addData(i, i, tmp * (1 - freep[1][num]));
                double rhoB = singlePipe.u[num] >= 0 ? singlePipe.rho[num] : singlePipe.rho[num + 1];
                double flux = rhoB * singlePipe.u[num] * area;
                sourcep[i] += flux;
                residualtmp += flux;
                double xn = singlePipe.dp[num];
                sourcep[i] += tmp * xn;
            }
        }
        csrp.addData(i, i, volume[i] / c2 / deltaT);
        double nonsteady = volume[i] * (rho0[i] - rho[i]) / deltaT;
        sourcep[i] += nonsteady;
        sourcep[i] += node.mass;
        residualtmp += nonsteady + node.mass;
        residual += residualtmp * residualtmp;
    }

    private void boundaryPressure(CSRMatrix csrp, double[] sourcep, int i, Node node) {
        ArrayList<Boolean> isSource = node.isSource;
        ArrayList<Integer> edgeID = node.edgeID;
        sourcep[i] = 0;
        int id = edgeID.get(0);
        Edge edge = edges[id];
        SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
        double[][] freep = edge.freep;
        int num = singlePipe.num;
        if (isSource.get(0)) {
            double tmp = singlePipe.cp[0];
            csrp.setData(i, edge.target, tmp * freep[1][1]);
            // tmp如果等于0，由说明为定压结点
            csrp.setData(i, i, tmp == 0 ? 1 : tmp * (freep[0][1] - 1));
            double x1 = singlePipe.dp[1];
            sourcep[i] += -tmp * x1;
        } else {
            double tmp = singlePipe.ap[num];
            csrp.setData(i, edge.source, -tmp * freep[0][num]);
            // tmp如果等于0，由说明为定压结点
            csrp.setData(i, i, tmp == 0 ? 1 : tmp * (1 - freep[1][num]));
            double xn = singlePipe.dp[num];
            sourcep[i] += tmp * xn;
        }
    }

    private void nodeMomentum(CSRMatrix csrv, double[] sourcev, int i, Node node) {
        ArrayList<Boolean> isSource = node.isSource;
        ArrayList<Integer> edgeID = node.edgeID;
        int n = edgeID.size();
        for (int j = 0; j < n; j++) {
            int id = edgeID.get(j);
            Edge edge = edges[id];
            SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
            double[][] freev = edge.freev;
            double deltaT = singlePipe.deltaT;
            if (isSource.get(j)) {// 如果此节点为此边的起点
                int row = 2 * id;
                sourcev[row] = singlePipe.sm[0];
                // 对流项采用迎风格式
                if (singlePipe.u[0] < 0) {// 管道内部
                    double b = singlePipe.bm[0];
                    double c = singlePipe.cm[0];
                    double x1 = singlePipe.u[1];
                    double f0 = freev[0][1];
                    double fn = freev[1][1];
                    csrv.addData(row, row, b + c * f0);
                    csrv.addData(row, row + 1, c * fn);
                    sourcev[row] += -c * x1;
                } else {// 节点处，需要使用节点离散方法
                    double rhoc = (singlePipe.rho[0] + singlePipe.rho[1]) / 2;
                    double b = rhoc / deltaT;
                    b += singlePipe.lambda[0] * rhoc * Math.abs(singlePipe.u[0]) / 2 / singlePipe.D;
                    csrv.addData(row, row, b);
                    double V = volume[i];
                    double ub = singlePipe.u[0];
                    double rhob = singlePipe.rho[0];
                    double tmp = ub * rhob;
                    // 节点对流项
                    for (int k = 0; k < n; k++) {
                        int idk = edgeID.get(k);
                        Edge edgek = edges[idk];
                        SinglePipeDiscretizer pipek = edgek.pipeDiscretizer;
                        boolean pos = isSource.get(k);
                        int direct = pos ? 1 : -1;// 外法线方向
                        int col = 2 * idk + (pos ? 0 : 1);// 节点位置
                        double Ak = pipek.Area;
                        csrv.addData(row, col, Ak * direct * tmp / V);
                    }
                }
            } else {// 如果此节点为此边的终点
                int row = 2 * id + 1;
                int num = singlePipe.num;
                sourcev[row] = singlePipe.sm[num];
                // 对流项采用迎风格式
                if (singlePipe.u[num] >= 0) {// 管道内部
                    double b = singlePipe.bm[num];
                    double a = singlePipe.am[num - 1];
                    double xn = singlePipe.u[num - 1];
                    double f0 = freev[0][num - 1];
                    double fn = freev[1][num - 1];
                    csrv.addData(row, row, b + a * fn);
                    csrv.addData(row, row - 1, a * f0);
                    sourcev[row] += -a * xn;
                } else {// 节点处，需要使用节点离散方法
                    double rhoc = (singlePipe.rho[num] + singlePipe.rho[num + 1]) / 2;
                    double b = rhoc / deltaT;
                    b += singlePipe.lambda[num] * rhoc * Math.abs(singlePipe.u[num]) / 2 / singlePipe.D;
                    csrv.addData(row, row, b);
                    double V = volume[i];
                    double ub = singlePipe.u[num];
                    double rhob = singlePipe.rho[num];
                    double tmp = ub * rhob;
                    // 节点对流项
                    for (int k = 0; k < n; k++) {
                        int idk = edgeID.get(k);
                        Edge edgek = edges[idk];
                        SinglePipeDiscretizer pipek = edgek.pipeDiscretizer;
                        boolean pos = isSource.get(k);
                        int direct = pos ? 1 : -1;// 外法线方向
                        int col = 2 * idk + (pos ? 0 : 1);// 节点位置
                        double Ak = pipek.Area;
                        csrv.addData(row, col, Ak * direct * tmp / V);
                    }
                }
            }
        }
    }

    private void boundaryMomentum(CSRMatrix csrv, double[] sourcev, Node node) {
        ArrayList<Boolean> isSource = node.isSource;
        ArrayList<Integer> edgeID = node.edgeID;
        int id = edgeID.get(0);
        Edge edge = edges[id];
        SinglePipeDiscretizer singlePipe = edge.pipeDiscretizer;
        double[][] freev = edge.freev;
        if (isSource.get(0)) {
            int row = 2 * id;
            double b = singlePipe.bm[0];
            double c = singlePipe.cm[0];
            double x1 = singlePipe.u[1];
            double f0 = freev[0][1];
            double fn = freev[1][1];
            csrv.setData(row, row, b + c * f0);
            csrv.setData(row, row + 1, c * fn);
            sourcev[row] = singlePipe.sm[0] - c * x1;
        } else {
            int row = 2 * id + 1;
            int num = singlePipe.num;
            double b = singlePipe.bm[num];
            double a = singlePipe.am[num - 1];
            double xn = singlePipe.u[num - 1];
            double f0 = freev[0][num - 1];
            double fn = freev[1][num - 1];
            csrv.addData(row, row, b + a * fn);
            csrv.addData(row, row - 1, a * f0);
            sourcev[row] = singlePipe.sm[num] - a * xn;
        }
    }

    private CSRMatrix createVelocityMatrix() {
        int numNode = nodes.length;
        int numEdge = edges.length;
        CSRMatrix.Builder csrBuilder = new CSRMatrix.Builder();
        csrBuilder.shape(numEdge * 2, numEdge * 2);
        for (int i = 0; i < numNode; i++) {
            Node node = nodes[i];
            ArrayList<Boolean> isSource = node.isSource;
            ArrayList<Integer> edgeID = node.edgeID;
            int n = edgeID.size();
            // 连接此节点的所有边上的点
            for (int j = 0; j < n; j++) {
                int id = edgeID.get(j);
                int row = 2 * id + (isSource.get(j) ? 0 : 1);
                // 连接此边节点的所有边节点，包括本身
                for (int k = 0; k < n; k++) {
                    int idk = edgeID.get(k);
                    int col = 2 * idk + (isSource.get(k) ? 0 : 1);
                    csrBuilder.addData(row, col, 0);
                }
                // 此边上的另一个节点
                int col = 2 * id + (isSource.get(j) ? 1 : 0);
                csrBuilder.addData(row, col, 0);
            }
        }
        return csrBuilder.build();
    }

    private CSRMatrix createPressureMatrix() {
        int numNode = nodes.length;
        CSRMatrix.Builder csrBuilder = new CSRMatrix.Builder();
        csrBuilder.shape(numNode, numNode);
        for (int i = 0; i < numNode; i++) {
            csrBuilder.addData(i, i, 0);
            Node node = nodes[i];
            ArrayList<Boolean> isSource = node.isSource;
            ArrayList<Integer> edgeID = node.edgeID;
            int n = edgeID.size();
            for (int j = 0; j < n; j++) {
                Edge edge = edges[edgeID.get(j)];
                if (isSource.get(j)) {
                    csrBuilder.addData(i, edge.target, 0);
                } else {
                    csrBuilder.addData(i, edge.source, 0);
                }
            }
        }
        return csrBuilder.build();
    }

    public static class NetworkBuilder {
        ArrayList<Node> nodes = new ArrayList<>();
        ArrayList<Edge> edges = new ArrayList<>();

        public NetworkBuilder addNode() {
            return addNode("Node" + nodes.size());
        }

        public NetworkBuilder addNode(String name) {
            Node node = new Node();
            node.name = name;
            nodes.add(node);
            return this;
        }

        public NetworkBuilder addEdge(int source, int target) {
            return addEdge("Edge" + edges.size(), source, target);
        }

        public NetworkBuilder addEdge(String name, int source, int target) {
            Edge edge = new Edge();
            int edgeid = edges.size();
            edge.name = name;
            edge.source = source;
            Node nodeS = nodes.get(source);
            nodeS.edgeID.add(edgeid);
            nodeS.isSource.add(true);
            edge.target = target;
            Node nodeT = nodes.get(target);
            nodeT.edgeID.add(edgeid);
            nodeT.isSource.add(false);
            edges.add(edge);
            return this;
        }

        public NetworkSolver build() {
            NetworkSolver ns = new NetworkSolver();
            int nds = nodes.size();
            int eds = edges.size();
            ns.volume = new double[nds];
            ns.dp = new double[nds];
            ns.rho = new double[nds];
            ns.T = new double[nds];
            ns.u = new double[eds * 2];
            ns.rho0 = new double[nds];
            ns.T0 = new double[nds];
            ns.nodes = nodes.toArray(new Node[0]);
            ns.edges = edges.toArray(new Edge[0]);
            return ns;
        }
    }

    public static class SolverBuilder {
        NetworkSolver solver;
        public SinglePipeSolver.Builder[] singleBuilders;
        int len;

        public SolverBuilder(NetworkSolver solver) {
            this.solver = solver;
            len = solver.edges.length;
            singleBuilders = new SinglePipeSolver.Builder[len];
            for (int i = 0; i < len; i++) {
                singleBuilders[i] = new SinglePipeSolver.Builder();
            }
        }

        public SolverBuilder fluidComponents(String[] names, double[] moleFractions) {
            for (int i = 0; i < len; i++) {
                singleBuilders[i].fluidComponents(names, moleFractions);
            }
            return this;
        }

        /**
         * @Author Shi Guoyun
         * @Description 管道环境设置，每个管道只能设置一个参数，可以通过获取单管进行详细设置
         * @Date 14:25:41 2021.01.19 019
         *
         * @param surroundingTemperature 环境温度，K
         * @param heatTransfers 总换热系数，W/(m2*K)
         */
        public SolverBuilder pipeSurrounding(double[] surroundingTemperature, double[] heatTransfers) {
            for (int i = 0; i < len; i++) {
                singleBuilders[i].pipeSurrounding(surroundingTemperature[i], null, heatTransfers[i], null);
            }
            return this;
        }

        public SolverBuilder discretize(double deltaT, double[] deltaX, int[] numOfCells) {
            for (int i = 0; i < len; i++) {
                singleBuilders[i].discretize(deltaT, deltaX[i], numOfCells[i]);
            }
            return this;
        }

        public SolverBuilder pipeConfig(double[] D, double[] roughness) {
            for (int i = 0; i < len; i++) {
                singleBuilders[i].pipeConfig(D[i], roughness[i]);
            }
            return this;
        }

        public SolverBuilder boundary(BoundaryCondition[] boundary) {
            for (int i = 0; i < len; i++) {
                singleBuilders[i].boundary(boundary[i]);
            }
            return this;
        }

        public SolverBuilder nodeMass(int index, double mass) {
            solver.nodes[index].mass = mass;
            return this;
        }

        public SolverBuilder nodePressure(int index, double pressure) {
            solver.nodes[index].pressure = pressure;
            return this;
        }

        public SolverBuilder nodeTemperature(int index, double temperature) {
            solver.nodes[index].temperature = temperature;
            return this;
        }

        public NetworkSolver build() {
            Edge[] edges = solver.edges;
            for (int i = 0; i < len; i++) {
                SinglePipeSolver singleSolver = singleBuilders[i].build();
                Edge edge = edges[i];
                edge.pipeDiscretizer = singleSolver.discretizer;
                edge.propertyUpdater = singleSolver.propertyUpdater;
                solver.deltaT = singleSolver.deltaT;
                int nump = singleSolver.discretizer.nump;
                int numv = singleSolver.discretizer.numv;
                edge.freep[0] = new double[nump];
                edge.freep[1] = new double[nump];
                edge.freet[0] = new double[nump];
                edge.freet[1] = new double[nump];
                edge.freev[0] = new double[numv];
                edge.freev[1] = new double[numv];
            }

            Node[] nodes = solver.nodes;
            int n = nodes.length;
            for (int i = 0; i < n; i++) {
                Node node = nodes[i];
                double V = node.edgeID.stream().map(j -> {
                    Edge edge = edges[j];
                    double a = edge.pipeDiscretizer.Area;
                    double dx = edge.pipeDiscretizer.deltaX;
                    return a * dx;
                }).reduce((a, b) -> a + b).get();
                solver.volume[i] = V / node.edgeID.size();
                // 节点定压
                if (node.pressure != 0) {
                    int m = node.edgeID.size();
                    for (int j = 0; j < m; j++) {
                        int id = node.edgeID.get(j);
                        Edge edge = edges[id];
                        double p = node.pressure;
                        if (node.isSource.get(j)) {
                            edge.pipeDiscretizer.p[0] = p;
                            double t = edge.pipeDiscretizer.T[0];
                            double den = edge.propertyUpdater.bwrs.getDensity(p / 1000, t);
                            double Mg = edge.propertyUpdater.bwrs.parameter.Mg;
                            den = den * Mg;
                            edge.pipeDiscretizer.rho[0] = den;
                        } else {
                            int end = edge.pipeDiscretizer.numv;
                            edge.pipeDiscretizer.p[end] = p;
                            double t = edge.pipeDiscretizer.T[end];
                            double den = edge.propertyUpdater.bwrs.getDensity(p / 1000, t);
                            double Mg = edge.propertyUpdater.bwrs.parameter.Mg;
                            den = den * Mg;
                            edge.pipeDiscretizer.rho[end] = den;
                        }
                    }
                } else {
                    int m = node.edgeID.size();
                    double p = 0, t = 0;
                    BWRS bwrs = edges[0].propertyUpdater.bwrs;
                    for (int j = 0; j < m; j++) {
                        int id = node.edgeID.get(j);
                        Edge edge = edges[id];
                        if (node.isSource.get(j)) {
                            p += edge.pipeDiscretizer.p[0];
                            t += edge.pipeDiscretizer.T[0];
                        } else {
                            int end = edge.pipeDiscretizer.numv;
                            p += edge.pipeDiscretizer.p[end];
                            t += edge.pipeDiscretizer.T[end];
                        }
                    }
                    p /= m;
                    t /= m;
                    double den = bwrs.getDensity(p / 1000, t);
                    den *= bwrs.parameter.Mg;
                    for (int j = 0; j < m; j++) {
                        int id = node.edgeID.get(j);
                        Edge edge = edges[id];
                        if (node.isSource.get(j)) {
                            edge.pipeDiscretizer.p[0] = p;
                            edge.pipeDiscretizer.p0[0] = p;
                            edge.pipeDiscretizer.T[0] = t;
                            edge.pipeDiscretizer.T0[0] = t;
                            edge.pipeDiscretizer.rho[0] = den;
                            edge.pipeDiscretizer.rho0[0] = den;
                        } else {
                            int end = edge.pipeDiscretizer.numv;
                            edge.pipeDiscretizer.p[end] = p;
                            edge.pipeDiscretizer.p0[end] = p;
                            edge.pipeDiscretizer.T[end] = t;
                            edge.pipeDiscretizer.T0[end] = t;
                            edge.pipeDiscretizer.rho[end] = den;
                            edge.pipeDiscretizer.rho0[end] = den;
                        }
                    }
                }
            }
            for (int i = 0; i < len; i++) {
                edges[i].propertyUpdater.propertyUpdate();
            }
            return solver;
        }
    }

    public static class Node {
        public String name;
        public double mass;// 质量流量kg/s
        public double pressure;// 压力，不为零时按定压处理
        public double temperature; // 温度，不为零时按定温处理
        public ArrayList<Integer> edgeID = new ArrayList<>();
        public ArrayList<Boolean> isSource = new ArrayList<>();
    }

    public static class Edge {
        public String name;
        public int source, target;
        public SinglePipeDiscretizer pipeDiscretizer;
        public PropertyUpdater propertyUpdater;
        public double[][] freep = new double[2][];
        public double[][] freev = new double[2][];
        public double[][] freet = new double[2][];
    }

}
