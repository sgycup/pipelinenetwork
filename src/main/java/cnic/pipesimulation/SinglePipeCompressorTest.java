package cnic.pipesimulation;

import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.27 027 20:41:32
 */
public class SinglePipeCompressorTest {

    public static void main(String[] args) throws IOException {
        String[] names = {"CH4", "C2H4", "C3H8", "N2", "CO2"};
        double[] fractions = {0.9707, 0.0017, 0.0002, 0.0071, 0.0203};
        BoundaryCondition boundary = new BoundaryCondition();
        boundary.isPressOut = true;
        boundary.isPressIn = true;
        boundary.isTempIn = true;
        boundary.pressure = 5e6;
        boundary.temperature = 30 + 273.15;

        SinglePipeSolver solver = new SinglePipeSolver.Builder()
                .fluidComponents(names, fractions)
                .discretize(200, 1000, 200)
                .pipeSurrounding(30 + 273.15, null, 0, null)
                .pipeConfig((590.0 / 2) * 1e-3, 0.025e-3)
                .boundary(boundary)
                .build();

        solver.observer = new Compressor(1);

        solver.simulate(86400 * 10);

        String result = Arrays.toString(solver.discretizer.p);
        result = result.substring(1, result.length() - 1);
        String result2 = Arrays.toString(solver.discretizer.u);
        result2 = result2.substring(1, result2.length() - 1);
        String result3 = Arrays.toString(solver.discretizer.T);
        result3 = result3.substring(1, result3.length() - 1);
        FileUtils.writeLines(new File("singlecompressortest.txt"), Arrays.asList(result, result2, result3));
    }
}
