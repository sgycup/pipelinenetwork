package simiulatonutil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * @Description:
 * @autor Shi Guoyun<sgy-421205643@163.com>
 * @date 2021.01.08 008 10:19:19
 */
public class CSRMatrix implements Matrix {

    public double[] data; // CSR format data array of the matrix
    public int[] indices; // CSR format index array of the matrix
    public int[] indptr; // CSR format index pointer array of the matrix
    public int rows, cols; // Shape of the matrix
    public int nnz; // Number of stored values, including explicit zeros

    @Override
    public void preprocess(double[] s){
        for (int i = 0; i < rows; i++) {
            double max = 0;
            for (int j = indptr[i]; j < indptr[i + 1]; j++) {
                double d = Math.abs(data[j]);
                if (max < d) {
                    max = d;
                }
            }
            for (int j = indptr[i]; j < indptr[i + 1]; j++) {
                data[j] /= max;
            }
            s[i] /= max;
        }
    }

    @Override
    public void multipy(double[] x, double[] y) {
        for (int i = 0; i < rows; i++) {
            double sum = 0;
            for (int j = indptr[i]; j < indptr[i + 1]; j++) {
                sum += data[j] * x[indices[j]];
            }
            y[i] = sum;
        }
    }

    public void setData(int row, int col, double value) {
        for (int j = indptr[row]; j < indptr[row + 1]; j++) {
            if (indices[j] == col) {
                data[j] = value;
                break;
            }
        }
    }

    public void addData(int row, int col, double value) {
        for (int j = indptr[row]; j < indptr[row + 1]; j++) {
            if (indices[j] == col) {
                data[j] += value;
                break;
            }
        }
    }

    public void clearData() {
        Arrays.fill(data, 0);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                boolean present = false;
                int index = 0;
                for (int k = indptr[i]; k < indptr[i + 1]; k++) {
                    if (j == indices[k]) {
                        index = k;
                        present = true;
                        break;
                    }
                }
                if (present) {
                    sb.append(data[index]).append("\t");
                } else {
                    sb.append("0.0\t");
                }
            }
            sb.append('\b').append("\n");
        }
        sb.append("]");
        return sb.toString();
    }

    public static class Builder {
        private int rows, cols;
        private int nnz = 0;
        private HashMap<Integer, ArrayList<Data>> data = new HashMap<>();
        public Builder shape(int rows, int cols) {
            this.rows = rows;
            this.cols = cols;
            return this;
        }

        /**
         * @param row 行号，从0开始
         * @param col 列号，从0开始
         * @param value 值
         * @return
         */
        public Builder addData(int row, int col, double value) {
            ArrayList<Data> d = data.computeIfAbsent(row, k -> new ArrayList<>());
            // 防止重复
            for (Data dd: d) {
                if (dd.index == col) {
                    dd.value += value;
                    return this;
                }
            }
            Data da = new Data();
            da.index = col;
            da.value = value;
            d.add(da);
            nnz++;
            return this;
        }
        public CSRMatrix build() {
            CSRMatrix csr = new CSRMatrix();
            csr.rows = rows;
            csr.cols = cols;
            csr.nnz = nnz;
            int[] indptr = new int[rows + 1];
            ArrayList<Integer> indices = new ArrayList<>();
            ArrayList<Double> values = new ArrayList<>();
            for (int i = 0; i < rows; i++) {
                ArrayList<Data> d = data.get(i);
                if (d != null) {
                    Data[] da = d.toArray(new Data[0]);
                    Arrays.sort(da);
                    for (Data dd : da) {
                        indices.add(dd.index);
                        values.add(dd.value);
                    }
                }
                indptr[i + 1] = values.size();
            }
            csr.indptr = indptr;
            csr.indices = new int[nnz];
            csr.data = new double[nnz];
            for (int i = 0; i < nnz; i++) {
                csr.indices[i] = indices.get(i);
                csr.data[i] = values.get(i);
            }
            return csr;
        }
        private static class Data implements Comparable {
            int index;
            double value;

            @Override
            public int compareTo(Object o) {
                Data other = (Data)o;
                return Integer.compare(index, other.index);
            }
        }
    }

    public static void main(String[] args) {
        CSRMatrix csr = new CSRMatrix.Builder()
                .shape(3, 3)
                .addData(1, 1, 1.0 / 3)
                .addData(1, 0, 2000)
                .addData(1, 0, 2000)
                .build();
        System.out.println(csr);
    }
}
