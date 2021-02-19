package cnic.pipesimulation;

import org.apache.commons.io.FileUtils;
import simiulatonutil.ParameterCalculator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.24 024 14:01:46
 */
public class FanDiEnergy {
    public static void main(String[] args) throws IOException {
        NetworkSolver.NetworkBuilder networkBuilder = new NetworkSolver.NetworkBuilder();
        for (int i = 0; i < 10; i++) {
            networkBuilder.addNode();
        }
        networkBuilder.addEdge(0, 1);
        networkBuilder.addEdge(1, 2);
        networkBuilder.addEdge(2, 3);
        networkBuilder.addEdge(3, 4);
        networkBuilder.addEdge(2, 5);
        networkBuilder.addEdge(5, 6);
        networkBuilder.addEdge(1, 7);
        networkBuilder.addEdge(3, 8);
        networkBuilder.addEdge(5, 9);
        NetworkSolver network = networkBuilder.build();
        int n = 9;
        double[] deltaX = new double[n];
        double[] diameter = new double[n];
        double[] roughness = new double[n];
        Arrays.fill(deltaX, 500);
        Arrays.fill(diameter, 0.988);
        Arrays.fill(roughness, 14e-3);
        int[] numOfCells = new int[]{83, 187, 15, 15, 60, 42, 1, 1, 1};
        for (int i = 0; i < n; i++) {
            numOfCells[i] *= 2;
        }
        String[] names = new String[]{"CH4", "C2H6", "C3H8", "i-C4H10", "n-C4H10", "i-C5H12",
                "n-C5H12", "C6H14", "N2", "CO2"};
        double[] fractions = new double[]{88.17, 5.532, 1.286, 0.224, 0.292, 0.105, 0.093,
                0.292, 2.052, 1.954};
        fractions = ParameterCalculator.weight2mole(names, fractions);
        double[] Ten = new double[n];
        double[] K = new double[n];
        Arrays.fill(Ten, 15 + 273.15);
        Arrays.fill(K, 1);
        BoundaryCondition[] boundary = new BoundaryCondition[n];
        for (int i = 0; i < n; i++) {
            boundary[i] = new BoundaryCondition();
            boundary[i].temperature = 24 + 273.15;
            boundary[i].pressure = 6.8e6;
        }
        boundary[0].pressure = 8.5e6;
        boundary[0].isPressIn = true;
        boundary[0].isTempIn = true;
        boundary[3].pressure = 7e6;
        boundary[3].isPressOut = true;
        boundary[5].pressure = 6.8e6;
        boundary[5].isPressOut = true;
        boundary[6].isFlowOut = true;
        boundary[6].mass = -4;
        boundary[7].isFlowOut = true;
        boundary[7].mass = -4;
        boundary[8].isFlowOut = true;
        boundary[8].mass = -30;
        NetworkSolver solver = new NetworkSolver.SolverBuilder(network)
                .fluidComponents(names, fractions)
                .discretize(200, deltaX, numOfCells)
                .pipeConfig(diameter, roughness)
                .pipeSurrounding(Ten, K)
                .boundary(boundary)
                .build();

        solver.EnergyStart = 86400 * 2;
        solver.simulate(86400 * 10);

        NetworkSolver.Edge[] edges = solver.edges;
        ArrayList<String> results = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            String result = Arrays.toString(edges[i].pipeDiscretizer.p);
            result = result.substring(1, result.length() - 1);
            results.add("pipe" + (i+1) + "\n" + result);
            result = Arrays.toString(edges[i].pipeDiscretizer.u);
            result = result.substring(1, result.length() - 1);
            results.add(result);
            result = Arrays.toString(edges[i].pipeDiscretizer.T);
            result = result.substring(1, result.length() - 1);
            results.add(result);
        }
        FileUtils.writeLines(new File("FanDi.txt"), results);
    }
}
