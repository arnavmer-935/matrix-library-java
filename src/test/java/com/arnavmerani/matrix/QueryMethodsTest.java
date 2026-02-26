package test.java.com.arnavmerani.matrix;

import main.java.com.arnavmerani.matrix.Matrix;
import main.java.com.arnavmerani.matrix.MatrixException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class QueryMethodsTest {

    private Matrix square3x3;
    private Matrix rectangular3x2;
    private Matrix identity3;
    private Matrix nearIdentity3;
    private Matrix notIdentity3;
    private Matrix upper3;
    private Matrix lower3;
    private Matrix notTriangular3;
    private Matrix zero3;
    private Matrix constant3;
    private Matrix symmetric3;
    private Matrix skewSymmetric3;
    private Matrix nonSymmetric3;
    private Matrix singular3;
    private Matrix nonSingular3;

    @BeforeEach
    void setup() {

        square3x3 = Matrix.ofRows(
                new double[]{2, 1, 3},
                new double[]{4, 5, 6},
                new double[]{7, 8, 10}
        );

        rectangular3x2 = Matrix.ofRows(
                new double[]{1, 2},
                new double[]{3, 4},
                new double[]{5, 6}
        );

        identity3 = Matrix.createIdentityMatrix(3);

        nearIdentity3 = Matrix.ofRows(
                new double[]{1.0000005, 0, 0},
                new double[]{0, 0.9999996, 0},
                new double[]{0, 0, 1.0000004}
        );

        notIdentity3 = Matrix.ofRows(
                new double[]{1, 0, 0},
                new double[]{0, 2, 0},
                new double[]{0, 0, 1}
        );

        upper3 = Matrix.ofRows(
                new double[]{5, 2, -1},
                new double[]{0, 3, 4},
                new double[]{0, 0, 6}
        );

        lower3 = Matrix.ofRows(
                new double[]{7, 0, 0},
                new double[]{2, 8, 0},
                new double[]{1, 3, 9}
        );

        notTriangular3 = Matrix.ofRows(
                new double[]{1, 2, 3},
                new double[]{4, 5, 6},
                new double[]{7, 8, 9}
        );

        zero3 = Matrix.ofRows(
                new double[]{0, 0, 0},
                new double[]{0, 0, 0},
                new double[]{0, 0, 0}
        );

        constant3 = Matrix.ofRows(
                new double[]{4, 4, 4},
                new double[]{4, 4, 4},
                new double[]{4, 4, 4}
        );

        symmetric3 = Matrix.ofRows(
                new double[]{2, -1, 3},
                new double[]{-1, 5, 0},
                new double[]{3, 0, 4}
        );

        skewSymmetric3 = Matrix.ofRows(
                new double[]{0, 2, -1},
                new double[]{-2, 0, 4},
                new double[]{1, -4, 0}
        );

        nonSymmetric3 = Matrix.ofRows(
                new double[]{1, 2, 3},
                new double[]{0, 5, 6},
                new double[]{7, 8, 9}
        );

        singular3 = Matrix.ofRows(
                new double[]{1, 2, 3},
                new double[]{2, 4, 6},
                new double[]{1, 1, 1}
        );

        nonSingular3 = Matrix.ofRows(
                new double[]{3, 1, 2},
                new double[]{0, 4, 5},
                new double[]{1, 0, 6}
        );
    }

    @Test
    @DisplayName("Square and rectangular detection works correctly.")
    void squareMatrixCheckWorks() {
        assertTrue(square3x3.isSquareMatrix());
        assertFalse(rectangular3x2.isSquareMatrix());
    }

    @Test
    @DisplayName("Identity detection respects tolerance.")
    void identityMatrixDetectionWorks() {
        assertTrue(identity3.isIdentityMatrix());
        assertTrue(nearIdentity3.isIdentityMatrix());
        assertFalse(notIdentity3.isIdentityMatrix());
        assertFalse(rectangular3x2.isIdentityMatrix());
    }

    @Test
    @DisplayName("Upper triangular detection works correctly.")
    void upperTriangularDetectionWorks() {
        assertTrue(upper3.isUpperTriangular());
        assertFalse(notTriangular3.isUpperTriangular());
        assertFalse(rectangular3x2.isUpperTriangular());
    }

    @Test
    @DisplayName("Lower triangular detection works correctly.")
    void lowerTriangularDetectionWorks() {
        assertTrue(lower3.isLowerTriangular());
        assertFalse(notTriangular3.isLowerTriangular());
        assertFalse(rectangular3x2.isLowerTriangular());
    }

    @Test
    @DisplayName("Zero and constant matrix detection works.")
    void zeroAndConstantChecksWork() {
        assertTrue(zero3.isZeroMatrix());
        assertTrue(constant3.isConstantMatrix());
        assertFalse(square3x3.isZeroMatrix());
        assertFalse(square3x3.isConstantMatrix());
    }

    @Test
    @DisplayName("Symmetric and skew-symmetric detection works.")
    void symmetryChecksWork() {
        assertTrue(symmetric3.isSymmetric());
        assertFalse(nonSymmetric3.isSymmetric());
        assertFalse(rectangular3x2.isSymmetric());

        assertTrue(skewSymmetric3.isSkewSymmetric());
        assertFalse(symmetric3.isSkewSymmetric());
        assertFalse(rectangular3x2.isSkewSymmetric());
    }

    @Test
    @DisplayName("Singular detection works correctly.")
    void singularCheckWorks() {
        assertTrue(singular3.isSingular());
        assertFalse(nonSingular3.isSingular());
        assertThrows(MatrixException.class, () -> rectangular3x2.isSingular());
    }
}