package cnic.pipesimulation;

import simiulatonutil.CompressorCurve;

import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.28 028 9:16:39
 */
public class Valve extends DefaultObserver {

    private int location;
    private double[] deltaP;
    private double[] p;
    private double[] u;
    private double[] am, bm, cm, sm;
    private double[] ap, cp;
    private double[] correctCoeff;
    private double[] at, bt, ct, st, T;
    private double Area, Rg;
    public double Cgmax = 200, Cgmin = 0;
    public double FR = 1;

    public Valve(int location) {
        this.location = location;
    }

    @Override
    public void simulationStart(SinglePipeSolver sovler) {
        this.deltaP = sovler.discretizer.deltaP;
        this.p = sovler.discretizer.p;
        this.u = sovler.discretizer.u;
        this.Area = sovler.discretizer.Area;
        this.am = sovler.discretizer.am;
        this.bm = sovler.discretizer.bm;
        this.cm = sovler.discretizer.cm;
        this.sm = sovler.discretizer.sm;
        this.ap = sovler.discretizer.ap;
        this.cp = sovler.discretizer.cp;
        this.at = sovler.discretizer.at;
        this.bt = sovler.discretizer.bt;
        this.ct = sovler.discretizer.ct;
        this.st = sovler.discretizer.st;
        this.T = sovler.discretizer.T;
        this.correctCoeff = sovler.discretizer.correctCoeff;
        this.Rg = sovler.propertyUpdater.bwrs.R * 1000;
        Rg = Rg / sovler.propertyUpdater.bwrs.parameter.Mg;
    }

    double time;
    @Override
    public void newTimestep(double time) {
        this.time = time;
        if (FR < 0.1) {
            return;
        }
        if (FR > 0.9) {
            deltaP[location] = 0;
            return;
        }
        double ui = u[location];
        double pi = ui >= 0 ? p[location] : p[location + 1];
        ui = Math.abs(ui);
        double Q = ui * Area * 3600;
        double Cg = Cgmin + FR * (Cgmax - Cgmin);
        double po = Q / Cg;
        po = po * po;
        pi /= 1000;
        po = po * pi / Rg / 1.209;
        po = Math.sqrt(pi * pi - po);
        if (Double.isNaN(po)) {
            po = 0.01 * pi;
        }
        double dp = po - pi;
        deltaP[location] = dp > 0 ? 0 : dp * 1000;
    }

    @Override
    public void afterDiscretizeMomentum() {
        if (FR < 0.1) {
            am[location - 1] = 0;
            bm[location] = 1;
            cm[location] = 0;
            sm[location] = 0;
        }
    }

    @Override
    public void afterDiscretizePressure() {
        if (FR < 0.1) {
            correctCoeff[location] = 0;
            cp[location] = 0;
            ap[location] = 0;
        }
    }

    @Override
    public void afterDiscretizeTemperature() {
        if (time >= 0) {
            Arrays.fill(at, 0);
            Arrays.fill(ct, 0);
            Arrays.fill(bt, 1);
            System.arraycopy(T, 0, st, 0, T.length);
        }
    }
}
