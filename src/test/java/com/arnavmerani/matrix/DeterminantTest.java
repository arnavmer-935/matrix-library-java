package com.arnavmerani.matrix;

import com.arnavmerani.matrix.Matrix;
import com.arnavmerani.matrix.MatrixException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class DeterminantTest {

    private Matrix m1x1;
    private Matrix m2x2;
    private Matrix identity4;
    private Matrix zeroRow3x3;
    private Matrix dependent3x3;
    private Matrix upperTriangular3x3;
    private Matrix swapA;
    private Matrix swapB;
    private Matrix nonSquare2x3;
    private Matrix smallFloating3x3;
    private Matrix ugly4x4;
    private Matrix multiSwap4x4;

    @BeforeEach
    void setup() {

        m1x1 = Matrix.ofRows(
                new double[]{5}
        );

        m2x2 = Matrix.ofRows(
                new double[]{3, 8},
                new double[]{4, 6}
        );

        identity4 = Matrix.ofRows(
                new double[]{1, 0, 0, 0},
                new double[]{0, 1, 0, 0},
                new double[]{0, 0, 1, 0},
                new double[]{0, 0, 0, 1}
        );

        zeroRow3x3 = Matrix.ofRows(
                new double[]{1, 2, 3},
                new double[]{0, 0, 0},
                new double[]{7, 8, 9}
        );

        dependent3x3 = Matrix.ofRows(
                new double[]{1, 2, 3},
                new double[]{2, 4, 6},
                new double[]{7, 8, 9}
        );

        upperTriangular3x3 = Matrix.ofRows(
                new double[]{2, 3, 1},
                new double[]{0, 5, 4},
                new double[]{0, 0, 7}
        );

        swapA = Matrix.ofRows(
                new double[]{1, 2},
                new double[]{3, 4}
        );

        swapB = Matrix.ofRows(
                new double[]{3, 4},
                new double[]{1, 2}
        );

        nonSquare2x3 = Matrix.ofRows(
                new double[]{1, 2, 3},
                new double[]{4, 5, 6}
        );

        smallFloating3x3 = Matrix.ofRows(
                new double[]{1, 1, 1},
                new double[]{1, 1.000001, 1},
                new double[]{1, 1, 1.000002}
        );

        ugly4x4 = Matrix.ofRows(
                new double[]{4.2, -3.1, 0.5, 2.0},
                new double[]{1.0, 6.3, -2.4, 3.3},
                new double[]{-7.2, 4.4, 5.5, -1.1},
                new double[]{2.2, 0.0, -3.3, 8.8}
        );

        multiSwap4x4 = Matrix.ofRows(
                new double[]{0, 2, 3, 1},
                new double[]{0, 0, 4, 2},
                new double[]{5, 6, 0, 3},
                new double[]{1, 2, 3, 4}
        );
    }

    @Test
    @DisplayName("Determinant of 1x1 matrix returns the element.")
    void singleElementDeterminant_ReturnsElement() {
        assertEquals(5.0, m1x1.determinant(), 1e-6);
    }

    @Test
    @DisplayName("Determinant of 2x2 matrix returns correct computation.")
    void determinantOfTwoByTwoIsCorrect() {
        assertEquals(-14.0, m2x2.determinant(), 1e-6);
    }

    @Test
    @DisplayName("Determinant of Identity Matrix is 1.")
    void identityMatrixDeterminantIsOne() {
        assertEquals(1.0, identity4.determinant(), 1e-6);
    }

    @Test
    @DisplayName("Determinant of a matrix with a zero row is 0.")
    void determinantOfMatrixWithZeroRowReturnsZero() {
        assertEquals(0.0, zeroRow3x3.determinant(), 1e-6);
    }

    @Test
    @DisplayName("Determinant of a matrix with linearly dependent rows is 0.")
    void determinantOfMatrixWithLinearlyDependentRowsReturnsZero() {
        assertEquals(0.0, dependent3x3.determinant(), 1e-6);
    }

    @Test
    @DisplayName("Determinant of a matrix in Row Echelon Form (REF) is the product of diagonal elements.")
    void upperTriangularMatrixUsesDiagonalProduct() {
        assertEquals(70.0, upperTriangular3x3.determinant(), 1e-6);
    }

    @Test
    @DisplayName("Swapping a row negates determinant.")
    void rowSwapFlipsDeterminantSign() {
        double detA = swapA.determinant();
        double detB = swapB.determinant();

        assertEquals(-detA, detB, 1e-6);
    }

    @Test
    @DisplayName("Throws for Non-square matrices.")
    void nonSquareMatrixThrowsException() {
        assertThrows(MatrixException.class, () -> nonSquare2x3.determinant());
    }

    @Test
    @DisplayName("Tolerance works for very small determinants.")
    void determinantOfSmallFloatingPointNumbers() {
        assertEquals(2e-12, smallFloating3x3.determinant(), 1e-6);
    }

    @Test
    @DisplayName("Testing consistency of determinant results.")
    void determinantOfUgly4x4() {
        double det1 = ugly4x4.determinant();
        double det2 = ugly4x4.determinant();

        assertEquals(1503.1808, det1, 1e-6);
        assertEquals(det1, det2, 1e-6);
    }

    @Test
    @DisplayName("Testing determinant operations for matrices involving multiple row swaps.")
    void determinantOfMultiSwap4x4() {
        double expected = 84.0;
        double actual = multiSwap4x4.determinant();

        assertEquals(expected, actual, 1e-6);
    }

}
