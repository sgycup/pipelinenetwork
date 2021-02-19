package cnic.pipesimulation;

import simiulatonutil.CompressorCurve;

import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.27 027 20:04:24
 */
public class Compressor extends DefaultObserver {

    private SinglePipeSolver sovler;
    private int location;
    private double[] deltaP;
    private double[] p;
    private double[] u;
    private double[] at, bt, ct, st, T;
    private double Area;
    private double A, B;

    public Compressor(int location) {
        this.location = location;
    }

    @Override
    public void simulationStart(SinglePipeSolver sovler) {
        this.sovler = sovler;
        this.deltaP = sovler.discretizer.deltaP;
        this.p = sovler.discretizer.p;
        this.u = sovler.discretizer.u;
        this.Area = sovler.discretizer.Area;
        this.at = sovler.discretizer.at;
        this.bt = sovler.discretizer.bt;
        this.ct = sovler.discretizer.ct;
        this.st = sovler.discretizer.st;
        this.T = sovler.discretizer.T;
        A = CompressorCurve.A[0];
        B = CompressorCurve.B[0];
    }

    double time;
    @Override
    public void newTimestep(double time) {
        this.time = time;
        double pi = p[location];
        double ui = u[location];
        double Q = ui * Area * 3600;
        double epsilon = Math.sqrt(A - B * Q * Q);
        double dp = pi * (epsilon - 1);
        deltaP[location] = dp < 0 ? 0 : dp;
    }

    @Override
    public void afterDiscretizeTemperature() {
        if (time < 86400 / 2) {
            Arrays.fill(at, 0);
            Arrays.fill(ct, 0);
            Arrays.fill(bt, 1);
            System.arraycopy(T, 0, st, 0, T.length);
        }
    }
}
