package com.marchofer.Gauss;

import java.util.Arrays;

public class Gauss {
    private static boolean debug = false;


    private static double[][] solve(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            if (isColZero(matrix, i)) continue;
            matrix = gauss(matrix, i);
            if (matrix == null) return null;
            print(matrix, "Gauss: " + i);
        }
        for (int i = 0; i < matrix.length; i++) {
            if (matrix[i][i] == 0.0) continue;
            double k = 1.0 / matrix[i][i];
            matrix = mult(matrix, i, k);
            print(matrix, "Norm: " + i);
        }
        for (int i = matrix.length - 1; i >= 0; i--) {
            if (isColZero(matrix, i)) continue;
            matrix = rref(matrix, i);
            if (matrix == null) return null;
            print(matrix, "Rref: " + i);
        }
        return matrix;
    }

    private static double[][] gauss(double[][] matrix, int col) {
        for (int i = col + 1; i < matrix.length; i++) {
            int b = col;
            if (matrix[b][col] == 0.0) {
                for (int j = 0; j < matrix.length; j++) {
                    if (j == i) continue;
                    b = j;
                    if (matrix[b][col] != 0.0) break;
                }
            }
            if (matrix[b][col] == 0.0) return null;
            double k = - matrix[i][col] / matrix[b][col];
            matrix = add(matrix, col, i, b, k);
        }
        return matrix;
    }

    private static double[][] rref(double[][] matrix, int col) {
        for (int i = col - 1; i >= 0; i--) {
            int b = col;
            if (matrix[b][col] == 0.0) {
                for (int j = 0; j < matrix.length; j++) {
                    if (j == i) continue;
                    b = j;
                    if (matrix[b][col] != 0.0) break;
                }
            }
            if (matrix[b][col] == 0.0) return null;
            double k = - matrix[i][col] / matrix[b][col];
            matrix = add(matrix, col, i, b, k);
        }
        return matrix;
    }

    private static double[][] add(double[][] matrix, int col, int i, int b, double k) {
        for (int j = 0; j < matrix[0].length; j++) {
            matrix[i][j] += k * matrix[b][j];
        }
        return matrix;
    }

    private static double[][] mult(double[][] matrix, int i, double k) {
        for (int j = 0; j < matrix[0].length; j++) {
            matrix[i][j] *= k;
        }
        return matrix;
    }

    private static boolean isColZero(double[][] matrix, int col) {
        boolean result = true;
        for (int i = 0; i < matrix.length; i++) {
            result &= (matrix[i][col] == 0.0);
        }
        return result;
    }

    private static boolean isRowZero(double[][] matrix, int row) {
        boolean result = true;
        for (int i = 0; i < matrix.length; i++) {
            result &= (matrix[row][i] == 0.0);
        }
        return result;
    }

    private static void print(double[][] matrix, String text) {
        if (!debug) return;
        System.out.println(text);
        for (int i = 0; i < matrix.length; i++) {
            System.out.println(Arrays.toString(matrix[i]));
        }
        System.out.println();
    }

    private static void printRes(double[][] matrix) {
        System.out.println("Result:");
        for (int i = 0; i < matrix.length; i++) {
            System.out.println("x_" + i + " = " + matrix[i][matrix[0].length - 1]);
        }
        System.out.println();
    }

    public static void main(String[] args) {
        double[][] matrix = new double[][] {
                {1.0, 1.0, 1.0, 1.0},
                {2.0, 3.0, 6.0, 2.0},
                {1.0, 5.0, 9.0, 3.0}
        };

        matrix = solve(matrix);
        if (matrix == null) System.out.println("Singulary Matrix");
        else printRes(matrix);

        if (matrix != null && isColZero(matrix, matrix[0].length - 1)) {
            System.out.println("Homogenous matrix");
        }
    }
}
