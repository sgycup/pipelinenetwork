package cnic.pipesimulation;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.27 027 19:40:10
 */
public interface SolverObserver {

    /**
     * @Author Shi Guoyun
     * @Description 模拟开始
     * @Date 19:49:24 2021.01.27 027
     * @return
     */
    void simulationStart(SinglePipeSolver sovler);

    /**
     * @Author Shi Guoyun
     * @Description 新时步模拟开始
     * @Date 19:50:47 2021.01.27 027
     * @return
     */
    void newTimestep(double time);

    /**
     * @Author Shi Guoyun
     * @Description 内循环的开始，物性更新前，每次内循环都会调用
     * @Date 19:51:19 2021.01.27 027
     * @return
     */
    void innerLoopStart();

    /**
     * @Author Shi Guoyun
     * @Description 动量方程离散前，物性更新后，每次内循环都会调用
     * @Date 19:51:59 2021.01.27 027
     * @return
     */
    void beforeDiscretizeMomentum();

    /**
     * @Author Shi Guoyun
     * @Description 动量方程离散后、求解前，每次内循环都会调用
     * @Date 19:52:51 2021.01.27 027
     * @return
     */
    void afterDiscretizeMomentum();

    /**
     * @Author Shi Guoyun
     * @Description 压力方程离散前，动量方程求解后，每次内循环都会调用
     * @Date 19:53:29 2021.01.27 027
     * @return
     */
    void beforeDiscretizePressure();

    /**
     * @Author Shi Guoyun
     * @Description 压力方程离散后、求解前，每次内循环都会调用
     * @Date 19:54:45 2021.01.27 027
     * @return
     */
    void afterDiscretizePressure();

    /**
     * @Author Shi Guoyun
     * @Description 变量修正前，压力方程求解后，每次内循环都会调用
     * @Date 19:55:16 2021.01.27 027
     * @return
     */
    void beforeCorrect();

    /**
     * @Author Shi Guoyun
     * @Description 温度方程离散前，变量修正后，每次内循环都会调用
     * @Date 19:56:50 2021.01.27 027
     * @return
     */
    void beforeDiscretizeTemperature();

    /**
     * @Author Shi Guoyun
     * @Description 温度方程离散后、求解前，每次内循环都会调用
     * @Date 19:57:21 2021.01.27 027
     * @return
     */
    void afterDiscretizeTemperature();

    /**
     * @Author Shi Guoyun
     * @Description 模拟结束
     * @Date 19:58:54 2021.01.27 027
     * @return
     */
    void simulationEnd();
}
