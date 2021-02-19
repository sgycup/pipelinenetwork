package cnic.pipesimulation;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.13 013 10:23:59
 */
public class BoundaryCondition {
    // 出入口定流标志
    public boolean isFlowOut, isFlowIn;
    // 出入口定压标志
    public boolean isPressOut, isPressIn;
    // 出入口定温标志
    public boolean isTempOut, isTempIn;
    // 流量(kg/s)，压力(Pa)，温度(K)。压力温度同时作为初始条件
    public double mass, pressure, temperature;
}
