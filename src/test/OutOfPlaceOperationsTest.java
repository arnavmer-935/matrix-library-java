package test;

import main.Matrix;
import main.MatrixException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class OutOfPlaceOperationsTest {
    Matrix A;
    Matrix B;
    Matrix C;
    Matrix rectangular, rectangular2;
    Matrix identity2x2;
    Matrix singular2x2;

    @BeforeEach
    void setup() {
        A = Matrix.ofRows(
                new double[]{1, 2},
                new double[]{3, 4}
        );

        B = Matrix.ofRows(
                new double[]{5, 6},
                new double[]{7, 8}
        );

        C = Matrix.ofRows(
                new double[]{2, -1},
                new double[]{0, 3}
        );

        rectangular = Matrix.ofRows(
                new double[]{1, 2, 3},
                new double[]{4, 5, 6}
        );

        rectangular2 = Matrix.ofRows(
                new double[]{6,5,4},
                new double[]{3,2,1}
        );

        identity2x2 = Matrix.createIdentityMatrix(2);

        singular2x2 = Matrix.ofRows(
                new double[]{1, 2},
                new double[]{2, 4}
        );
    }

    @Nested
    class ScalarMultiplication {

        @Test
        @DisplayName("Out of place scalar multiplication occurs correctly.")
        void doesScalarOutOfPlaceWork() {
            Matrix result1 = A.multiplyByScalar(2.0);
            Matrix expected = new Matrix(new double[][]{
                                            {2,4},
                                            {6,8}}
            );
            assertEquals(result1, expected);

            Matrix result2 = A.multiplyByScalar(0.0);
            Matrix expected2 = Matrix.zeroMatrix(2, 2);
            assertEquals(result2, expected2);

            assertThrows(MatrixException.class, () -> A.multiplyByScalar(Double.NaN));
            assertThrows(MatrixException.class, () -> A.multiplyByScalar(Double.POSITIVE_INFINITY));
            assertThrows(MatrixException.class, () -> A.multiplyByScalar(Double.NEGATIVE_INFINITY));

            Matrix result3 = A.multiplyByScalar(1.0);
            assertEquals(result3, A);

            Matrix result4 = A.multiplyByScalar(-2);
            Matrix expected4 = new Matrix(new double[][] {{-2,-4},{-6,-8}});
            assertEquals(result4, expected4);
        }

        @Test
        @DisplayName("Out of place scalar multiplication does not mutate original operand.")
        void scalarMultiply_doesNotMutateOriginal() {
            Matrix copy = new Matrix(A);
            Matrix newCopy = A.multiplyByScalar(3.0);
            assertEquals(copy, A);
        }

        @Test
        @DisplayName("Out of place scalar multiplication returns a new Matrix instance.")
        void scalarMultiply_returnsNewInstance() {
            Matrix result = A.multiplyByScalar(2.0);
            assertNotSame(A, result);
        }
    }

    @Nested
    class OutOfPlaceAddition {

        @Test
        @DisplayName("Tests ideal path for out of place addition.")
        void doesAddOutOfPlaceWork() {

            Matrix idealPathResult = A.add(B);
            Matrix expected = Matrix.ofRows(new double[] {6,8}, new double[] {10,12});
            assertEquals(expected, idealPathResult);

            Matrix addingToNegativeSelf = A.add(A.multiplyByScalar(-1));
            Matrix expected2 = Matrix.zeroMatrix(2,2);
            assertEquals(expected2, addingToNegativeSelf);

            Matrix addingToZeroMatrix = A.add(Matrix.zeroMatrix(2,2));
            Matrix expected3 = new Matrix(A);
            assertEquals(expected3, addingToZeroMatrix);

            Matrix addingToSelf = A.add(A);
            Matrix expected4 = Matrix.ofRows(new double[]{2,4}, new double[]{6,8});
            assertEquals(expected4, addingToSelf);

            Matrix addRectangular = rectangular.add(rectangular2);
            Matrix expected5 = Matrix.constant(2,3,7);
            assertEquals(expected5, addRectangular);
        }

        @Test
        @DisplayName("Tests handling of invalid input.")
        void doesAddThrowForBadInput() {

            assertThrows(MatrixException.class, () -> {
                Matrix diffDimensions = A.add(rectangular);
            });
        }
    }
}
