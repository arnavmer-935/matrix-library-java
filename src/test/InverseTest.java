package test;

import main.Matrix;
import main.MatrixException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class InverseTest {

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

        identity4 = Matrix.createIdentityMatrix(4);

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
    @DisplayName("Inverse of 1x1 matrix returns reciprocal.")
    void singleElementInverse_ReturnsReciprocal() {
        Matrix inverse = m1x1.inverse();
        assertEquals(0.2, inverse.getEntry(0, 0), 1e-6);
    }

    @Test
    @DisplayName("2x2 matrix inverse is correct.")
    void inverseOfTwoByTwoIsCorrect() {
        Matrix inv = m2x2.inverse();
        Matrix product = m2x2.multiply(inv);

        assertTrue(product.isIdentityMatrix());
    }

    @Test
    @DisplayName("Inverse of Identity Matrix is Identity.")
    void identityMatrixInverseIsIdentity() {
        Matrix inv = identity4.inverse();
        assertEquals(inv, identity4);
    }

    @Test
    @DisplayName("Matrix with zero row throws exception.")
    void inverseOfMatrixWithZeroRowThrows() {
        assertThrows(MatrixException.class, () -> zeroRow3x3.inverse());
    }

    @Test
    @DisplayName("Matrix with linearly dependent rows throws exception.")
    void inverseOfDependentMatrixThrows() {
        assertThrows(MatrixException.class, () -> dependent3x3.inverse());
    }

    @Test
    @DisplayName("Upper triangular matrix inverse satisfies A * A^-1 = I.")
    void upperTriangularMatrixInverseWorks() {
        Matrix inv = upperTriangular3x3.inverse();
        Matrix product = upperTriangular3x3.multiply(inv);

        assertTrue(product.isIdentityMatrix());
    }

    @Test
    @DisplayName("Row swaps do not affect invertibility correctness.")
    void swappedMatrixStillInvertsCorrectly() {
        Matrix invA = swapA.inverse();
        Matrix invB = swapB.inverse();

        assertTrue(swapA.multiply(invA).isIdentityMatrix());
        assertTrue(swapB.multiply(invB).isIdentityMatrix());
    }

    @Test
    @DisplayName("Non-square matrices throw exception.")
    void nonSquareMatrixThrowsException() {
        assertThrows(MatrixException.class, () -> nonSquare2x3.inverse());
    }

    @Test
    @DisplayName("Small floating point matrix inverse is numerically stable.")
    void inverseHandlesFloatingPointStability() {
        Matrix inv = smallFloating3x3.inverse();
        Matrix product = smallFloating3x3.multiply(inv);

        assertTrue(product.isIdentityMatrix());
    }

    @Test
    @DisplayName("Inverse consistency check on ugly 4x4.")
    void inverseOfUgly4x4IsConsistent() {
        Matrix inv1 = ugly4x4.inverse();
        Matrix inv2 = ugly4x4.inverse();

        assertEquals(inv1, inv2);
        assertTrue(ugly4x4.multiply(inv1).isIdentityMatrix());
    }

    @Test
    @DisplayName("Inverse works for matrix requiring multiple row swaps.")
    void inverseOfMultiSwap4x4Works() {
        Matrix inv = multiSwap4x4.inverse();
        Matrix product = multiSwap4x4.multiply(inv);

        assertTrue(product.isIdentityMatrix());
    }
}
