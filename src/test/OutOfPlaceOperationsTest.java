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
    Matrix D;
    Matrix E;
    Matrix E2;
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

        D = Matrix.ofRows(
                new double[] {1,2,3},
                new double[] {4,5,6},
                new double[] {7,8,9}
        );

        E = new Matrix(new double[][] {{1,2}});

        E2 = Matrix.ofRows(new double[] {1});

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
            Matrix expected2 = new Matrix(2, 2);
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
            Matrix expected2 = new Matrix(2,2);
            assertEquals(expected2, addingToNegativeSelf);

            Matrix addingToZeroMatrix = A.add(new Matrix(2,2));
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

    @Nested
    class OutOfPlaceSubtraction {
        @Test
        void subtract_valid() {
            Matrix result = B.subtract(A);

            Matrix expected = Matrix.ofRows(
                    new double[]{4, 4},
                    new double[]{4, 4}
            );

            assertEquals(expected, result);
        }

        @Test
        void subtract_self_returnsZero() {
            Matrix result = A.subtract(A);
            Matrix expected = new Matrix(2,2);
            assertEquals(expected, result);
        }

        @Test
        void subtract_throws() {
            assertThrows(MatrixException.class, () -> A.subtract(rectangular));
            assertThrows(IllegalArgumentException.class, () -> A.subtract(null));
        }

        @Test
        void subtract_doesNotMutateOperands() {
            Matrix copyA = new Matrix(A);
            Matrix copyB = new Matrix(B);

            Matrix result = A.subtract(B);

            assertEquals(copyA, A);
            assertEquals(copyB, B);
        }
    }

    @Nested
    class OutOfPlaceTranspose {

        @Test
        @DisplayName("Transposes 2x2 and 3x3 square matrices and returns a new instance.")
        void transposeWorksFor2x2and3x3_andReturnsNewInstance() {

            Matrix transposeA = A.transpose();
            Matrix transposeD = D.transpose();

            Matrix expected1 = Matrix.ofRows(
                    new double[]{1,3},
                    new double[] {2,4}
            );

            Matrix expected2 = Matrix.ofRows(
                    new double[] {1,4,7},
                    new double[] {2,5,8},
                    new double[] {3,6,9}
            );

            assertEquals(expected1, transposeA);
            assertEquals(expected2, transposeD);

            assertNotSame(transposeA, A); //checks that new object is created, so references should differ
            assertNotSame(transposeD, D);
        }

        @Test
        @DisplayName("Rectangular matrix transposes correctly.")
        void transposeWorksForRectangularMatrix() {
            Matrix rectangularTranspose1 = rectangular.transpose();
            Matrix rectangularTranspose2 = rectangular2.transpose();
            Matrix edgeCaseTranspose = E.transpose();
            Matrix edgeCaseTranspose2 = E2.transpose();


            Matrix expected1 = Matrix.ofRows(
                    new double[] {1,4},
                    new double[] {2,5},
                    new double[] {3,6}
            );

            Matrix expected2 = Matrix.ofRows(
                    new double[] {6,3},
                    new double[] {5,2},
                    new double[] {4,1}
            );

            Matrix expected3 = Matrix.ofRows(
                    new double[] {1},
                    new double[] {2}
            );

            Matrix expected4 = new Matrix(E2); //singleton matrix should remain unchanged

            assertEquals(expected1, rectangularTranspose1);
            assertEquals(expected2, rectangularTranspose2);
            assertEquals(expected3, edgeCaseTranspose);
            assertEquals(expected4, edgeCaseTranspose2);
        }

        @Test
        @DisplayName("Transpose works for special matrix types.")
        void transposeWorksForSpecialMatrices() {
            Matrix specialTranspose1 = Matrix.createIdentityMatrix(3).transpose();
            Matrix specialTranspose2 = Matrix.createScalarMatrix(3, 5).transpose();
            Matrix specialTranspose3 = new Matrix(3,3).transpose();
            Matrix specialTranspose4 = Matrix.constant(2,3,5).transpose();

            Matrix expected1 = Matrix.createIdentityMatrix(3);
            Matrix expected2 = Matrix.createScalarMatrix(3,5);
            Matrix expected3 = new Matrix(3,3);
            Matrix expected4 = Matrix.constant(3,2,5);

            assertEquals(expected1, specialTranspose1);
            assertEquals(expected2, specialTranspose2);
            assertEquals(expected3, specialTranspose3);
            assertEquals(expected4, specialTranspose4);
        }

        @Test
        @DisplayName("Immutability of original matrix during transposition.")
        void transposeLeavesOriginalUnchanged() {
            Matrix originalA = new Matrix(A);

            Matrix transposedA = originalA.transpose();

            Matrix expectedTranspose = Matrix.ofRows(
                    new double[] {1, 3},
                    new double[] {2, 4}
            );

            Matrix expectedOriginal = Matrix.ofRows(
                    new double[] {1, 2},
                    new double[] {3, 4}
            );

            assertNotSame(originalA, transposedA);
            assertEquals(expectedTranspose, transposedA);
            assertEquals(expectedOriginal, originalA);

            transposedA.setEntry(999, 0, 0);
            assertEquals(1, originalA.getEntry(0, 0));
        }

    }

    @Nested
    class MultiplicationTest {

        @Test
        @DisplayName("Multiplication works correctly for 2x2 matrices")
        void multiplicationTestFor2x2() {

            Matrix result1 = A.multiply(B);
            Matrix expected1 = Matrix.ofRows(
                    new double[]{19, 22},
                    new double[]{43, 50}
            );

            assertEquals(2, result1.getRows());
            assertEquals(2, result1.getColumns());
            assertEquals(expected1, result1);

            Matrix result2 = A.multiply(A);
            Matrix expected2 = Matrix.ofRows(
                    new double[]{7, 10},
                    new double[]{15, 22}
            );

            assertEquals(expected2, result2);
        }

        @Test
        @DisplayName("Multiplication works correctly for 3x3 matrices")
        void multiplicationTestFor3x3() {

            Matrix result = D.multiply(D);

            Matrix expected = Matrix.ofRows(
                    new double[]{30, 36, 42},
                    new double[]{66, 81, 96},
                    new double[]{102, 126, 150}
            );

            assertEquals(3, result.getRows());
            assertEquals(3, result.getColumns());
            assertEquals(expected, result);
        }

        @Test
        @DisplayName("Rectangular multiplication produces correct dimensions and values")
        void rectangularValidMultiplication() {

            Matrix result = rectangular.multiply(rectangular2.transpose());

            Matrix expected = Matrix.ofRows(
                    new double[]{28, 10},
                    new double[]{73, 28}
            );

            assertEquals(2, result.getRows());
            assertEquals(2, result.getColumns());
            assertEquals(expected, result);
        }

        @Test
        @DisplayName("Identity matrix does not change matrix when multiplied")
        void identityMultiplication() {

            Matrix left = identity2x2.multiply(A);
            Matrix right = A.multiply(identity2x2);

            assertEquals(A, left);
            assertEquals(A, right);
        }

        @Test
        @DisplayName("Multiplication by zero matrix produces zero matrix")
        void zeroMatrixMultiplication() {

            Matrix zero = Matrix.ofRows(
                    new double[]{0, 0},
                    new double[]{0, 0}
            );

            Matrix result = A.multiply(zero);

            Matrix expected = Matrix.ofRows(
                    new double[]{0, 0},
                    new double[]{0, 0}
            );

            assertEquals(expected, result);
        }

        @Test
        @DisplayName("Matrix multiplication is not commutative")
        void nonCommutativityTest() {

            Matrix ab = A.multiply(C);
            Matrix ba = C.multiply(A);

            assertNotEquals(ab, ba);
        }

        @Test
        @DisplayName("Invalid dimension multiplication throws exception")
        void invalidDimensionMultiplication() {

            assertThrows(MatrixException.class, () -> A.multiply(D));

            assertThrows(MatrixException.class,
                    () -> rectangular.multiply(rectangular2));
        }

        @Test
        @DisplayName("Multiplying with null throws exception")
        void nullMultiplicationThrows() {

            assertThrows(IllegalArgumentException.class,
                    () -> A.multiply(null));
        }

        @Test
        @DisplayName("Out-of-place multiplication does not mutate operands")
        void immutabilityTest() {

            Matrix originalA = Matrix.ofRows(
                    new double[]{1, 2},
                    new double[]{3, 4}
            );

            A.multiply(B);

            assertEquals(originalA, A);
            assertEquals(Matrix.ofRows(
                    new double[]{5, 6},
                    new double[]{7, 8}
            ), B);
        }
    }

}
