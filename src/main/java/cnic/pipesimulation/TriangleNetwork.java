package cnic.pipesimulation;

import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.14 014 13:03:32
 */
public class TriangleNetwork {
    public static void main(String[] args) throws IOException {
        NetworkSolver.NetworkBuilder networkBuilder = new NetworkSolver.NetworkBuilder();
        NetworkSolver network = networkBuilder
                .addNode()
                .addNode()
                .addNode()
                .addEdge(0, 1)
                .addEdge(1, 2)
                .addEdge(2, 0)
                .build();

        String[] names = {"CH4", "C2H4", "C3H8", "N2", "CO2"};
        double[] fractions = {0.9707, 0.0017, 0.0002, 0.0071, 0.0203};

        NetworkSolver.SolverBuilder solverBuilder = new NetworkSolver.SolverBuilder(network);
        solverBuilder.fluidComponents(names, fractions);

        double den = solverBuilder.singleBuilders[0].bwrs.getDensity(101.325, 273.15);
        double Mg = solverBuilder.singleBuilders[0].bwrs.parameter.Mg;
        den = den * Mg;
        BoundaryCondition boundary1 = new BoundaryCondition();
        boundary1.isPressIn = true;
        boundary1.pressure = 5e6;
        boundary1.temperature = 30 + 273.15;
        BoundaryCondition boundary2 = new BoundaryCondition();
        boundary2.pressure = 5e6;
        boundary2.temperature = 30 + 273.15;
        BoundaryCondition boundary3 = new BoundaryCondition();
        boundary3.isPressOut = true;
        boundary3.pressure = 5e6;
        boundary3.temperature = 30 + 273.15;

        NetworkSolver solver = solverBuilder
                .discretize(10, new double[]{500, 500, 500}, new int[]{180, 200, 160})
                .pipeConfig(new double[]{0.6, 0.6, 0.6}, new double[]{0.025e-3, 0.025e-3, 0.025e-3})
                .boundary(new BoundaryCondition[]{boundary1, boundary2, boundary3})
                .nodeMass(1, -15 * den)
                .nodeMass(2, -50 * den)
                .nodePressure(0, 5e6)
                .nodeTemperature(0, 30 + 273.15)
                .build();

        solver.simulate(3600 * 24);

        String result = Arrays.toString(solver.edges[0].pipeDiscretizer.p);
        result = result.substring(1, result.length() - 1);
        String result2 = Arrays.toString(solver.edges[0].pipeDiscretizer.u);
        result2 = result2.substring(1, result2.length() - 1);
        String result3 = Arrays.toString(solver.edges[0].pipeDiscretizer.T);
        result3 = result3.substring(1, result3.length() - 1);
        FileUtils.writeLines(new File("triangle.txt"), Arrays.asList(result, result2, result3));
    }
}
