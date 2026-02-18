package test;

import main.Matrix;

import main.MatrixException;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class InPlaceOperationsTest {

    @Nested
    class InPlaceScalar {
        @Test
        @DisplayName("Testing in-place scalar multiplication with zero.")
        void multiplyByScalarInPlace_zero() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2},
                    new double[]{3, 4}
            );

            A.multiplyByScalarInPlace(0.0);

            Matrix expected = new Matrix(2,2);

            assertEquals(expected, A);
        }

        @Test
        @DisplayName("Testing in-place scalar multiplication with 1.")
        void multiplyByScalarInPlace_one() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2},
                    new double[]{3, 4}
            );

            Matrix B = new Matrix(A);

            A.multiplyByScalarInPlace(1.0);
            assertEquals(B, A);
        }

        @Test
        @DisplayName("Testing in-place scalar multiplication with nonzero scalar.")
        void multiplyByScalarInPlace_regularScalar() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2},
                    new double[]{3, 4}
            );

            A.multiplyByScalarInPlace(2.0);

            Matrix expected = Matrix.ofRows(
                    new double[]{2, 4},
                    new double[]{6, 8}
            );

            assertEquals(expected, A);
        }

        @Test
        @DisplayName("Testing in-place scalar multiplication with NaN and infinite values.")
        void multiplyByScalarInPlace_undefinedScalar() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2},
                    new double[]{3, 4}
            );

            assertThrows(MatrixException.class, () -> A.multiplyByScalarInPlace(Double.NaN));
            assertThrows(MatrixException.class, () -> A.multiplyByScalarInPlace(Double.POSITIVE_INFINITY));
            assertThrows(MatrixException.class, () -> A.multiplyByScalarInPlace(Double.NEGATIVE_INFINITY));
        }
    }

    @Nested
    class InPlaceAddition {

        @Test
        @DisplayName("Testing ideal case of in-place addition.")
        void addInPlace_regularCase() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2},
                    new double[]{3, 4}
            );

            Matrix B = Matrix.ofRows(
                    new double[]{5, 6},
                    new double[]{7, 8}
            );

            A.addInPlace(B);

            Matrix expected = Matrix.ofRows(
                    new double[]{6, 8},
                    new double[]{10, 12}
            );

            assertEquals(expected, A);
        }

        @Test
        @DisplayName("Testing null case of in-place addition.")
        void addInPlace_nullThrows() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2}
            );

            assertThrows(IllegalArgumentException.class, () -> A.addInPlace(null));
        }

        @Test
        @DisplayName("Testing dimension mismatch case of in-place addition.")
        void addInPlace_orderMismatchThrows() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2}
            );

            Matrix B = Matrix.ofRows(
                    new double[]{1},
                    new double[]{2}
            );

            assertThrows(MatrixException.class, () -> A.addInPlace(B));
        }

        @Test
        @DisplayName("Testing mutation only of 'this' and not 'other' during in-place addition.")
        void addInPlace_doesNotMutateOther() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 1}
            );

            Matrix B = Matrix.ofRows(
                    new double[]{2, 2}
            );

            Matrix copyOfB = new Matrix(B);

            A.addInPlace(B);

            assertEquals(copyOfB, B);
        }
    }

    @Nested
    class InPlaceSubtraction {

        @Test
        @DisplayName("Tests ideal case of in-place subtraction.")
        void subtractInPlace_regularCase() {
            Matrix A = Matrix.ofRows(
                    new double[]{5, 6},
                    new double[]{7, 8}
            );

            Matrix B = Matrix.ofRows(
                    new double[]{1, 2},
                    new double[]{3, 4}
            );

            A.subtractInPlace(B);

            Matrix expected = Matrix.ofRows(
                    new double[]{4, 4},
                    new double[]{4, 4}
            );

            assertEquals(expected, A);
        }

        @Test
        @DisplayName("Test null argument case of in-place subtraction.")
        void subtractInPlace_nullThrows() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2}
            );

            assertThrows(IllegalArgumentException.class, () -> A.subtractInPlace(null));
        }

        @Test
        @DisplayName("Tests dimension mismatch case of in-place subtraction.")
        void subtractInPlace_orderMismatchThrows() {
            Matrix A = Matrix.ofRows(
                    new double[]{1, 2}
            );

            Matrix B = Matrix.ofRows(
                    new double[]{1},
                    new double[]{2}
            );

            assertThrows(MatrixException.class, () -> A.subtractInPlace(B));
        }

        @Test
        @DisplayName("Tests mutation of only 'this' and not 'other' in-place subtraction.")
        void subtractInPlace_doesNotMutateOther() {
            Matrix A = Matrix.ofRows(
                    new double[]{5, 5}
            );

            Matrix B = Matrix.ofRows(
                    new double[]{2, 3}
            );

            Matrix copyOfB = new Matrix(B);

            A.subtractInPlace(B);

            assertEquals(copyOfB, B);
        }
    }

    @Nested
    class TransposeInPlace {

        @Test
        @DisplayName("Transposes square matrices correctly.")
        void doesInPlaceTransposeWork() {
            Matrix A = Matrix.ofRows(
                    new double[] {-1,-2},
                    new double[] {-3,-4}
            );
            A.transposeInPlace();
            Matrix expected1 = Matrix.ofRows(
                    new double[] {-1,-3},
                    new double[] {-2,-4}
            );
            assertEquals(expected1, A);

            Matrix I3x3 = Matrix.createIdentityMatrix(3);
            I3x3.transposeInPlace();
            Matrix expected2 = Matrix.ofRows(
                    new double[] {1,0,0},
                    new double[] {0,1,0},
                    new double[] {0,0,1}
            );
            assertEquals(expected2, I3x3);

            Matrix B = Matrix.ofRows(
                    new double[] {1,2,3},
                    new double[] {4,5,6},
                    new double[] {7,8,9}
            );
            B.transposeInPlace();
            Matrix expected3 = Matrix.ofRows(
                    new double[]{1, 4, 7},
                    new double[]{2, 5, 8},
                    new double[]{3, 6, 9}
            );
            assertEquals(expected3, B);

            Matrix singleton = new Matrix(new double[][]{{1}});
            singleton.transposeInPlace();
            Matrix expected4 = Matrix.ofRows(new double[]{1});
            assertEquals(expected4, singleton);

        }

        @Test
        @DisplayName("In place matrix transposition throws for rectangular matrix")
        void transposeInPlaceThrowsForRectangular() {
            Matrix rectangular = Matrix.ofRows(
                    new double[]{1, 2, 3},
                    new double[]{4, 5, 6}
            );

            assertThrows(MatrixException.class, rectangular::transposeInPlace);
        }
    }


}
