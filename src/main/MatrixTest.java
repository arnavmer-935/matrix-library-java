package main;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

/*
PHASE 1 OF UNIT TESTING FOR THE FUNCTIONALITIES OF THE MATRIX CLASS.
 */
public class MatrixTest {

    @Nested
    @DisplayName("Creation & Factory Methods")
    class CreationAndFactoryTest {

        @Test
        @DisplayName("Row/column constructor creates zero-filled matrix with correct dimensions")
        void constructorWithRowsAndColsCreatesZeroMatrix() {

        }

        @Test
        @DisplayName("Square constructor creates matrix with equal row and column count")
        void squareConstructorCreatesCorrectDimensions() {

        }

        @Test
        @DisplayName("Constructor rejects negative or zero dimensions")
        void constructorRejectsInvalidDimensions() {

        }

        @Test
        @DisplayName("Constructor from 2D array performs deep copy")
        void constructorFrom2DArrayPerformsDeepCopy() {

        }

        @Test
        @DisplayName("Constructor rejects jagged 2D arrays")
        void constructorRejectsJaggedArrays() {

        }

        @Test
        @DisplayName("Constant factory fills matrix with the given value")
        void constantFactoryFillsMatrixWithValue() {

        }

        @Test
        @DisplayName("Scalar factory creates 1x1 matrix with the given value")
        void scalarFactoryCreatesUnitMatrix() {

        }

        @Test
        @DisplayName("ofRows method creates matrix correctly.")
        void ofRowsFactoryMethodCreatesCorrectMatrix() {

        }

        @Test
        @DisplayName("ofColumns method creates matrix correctly.")
        void ofColumnsFactoryMethodCreatesCorrectMatrix() {

        }
    }


    @Nested
    @DisplayName("Matrix Accessors and Mutators.")
    class GetterSetterTest {

        @Test
        @DisplayName("getOrder returns a correctly formatted tuple of matrix dimensions.")
        void getOrderWorks() {

        }

        @Test
        @DisplayName("setEntry method throws a MatrixException when trying to update a matrix entry.")
        void setEntryThrowsExceptionForOutOfBounds() {

        }

        @Test
        @DisplayName("setRow works")
        void setRowUpdatesMatrixCorrectly() {

        }

        @Test
        @DisplayName("setCol works")
        void setColUpdatesMatrixCorrectly() {

        }


    }

    @Nested
    @DisplayName("In-place Matrix Methods.")
    class InPlaceOperationsTest {
        //TODO: Stub this
    }

    @Nested
    @DisplayName("Out of place Matrix Methods.")
    class OutOfPlaceOperationsTest {
        //TODO: Stub this
    }

    @Nested
    @DisplayName("Numerical Matrix Methods.")
    class NumericalMethodsTest {
        //TODO: Stub this
    }

    @Nested
    @DisplayName("Matrix Query Methods.")
    class QueryMethodsTest {

        @Test
        @DisplayName("Checks whether a matrix is square.")
        void doesSquareWorkForTinyAndLargeMatrices() {

        }

        @Test
        @DisplayName("Checks whether isUpperTriangular works for a lower triangle with close-to-zero values.")
        void doesUpperTriangularWorkForNearlyZeroLT() {

        }

        @Test
        @DisplayName("Converse of isLowerTriangular")
        void doesLowerTriangularWorkForNearlyZeroUT() {

        }

        @Test
        @DisplayName("Checks whether isConstant works for entries close to one of the values.")
        void doesIsConstantMatrixWorkForCloseValues() {

        }

        @Test
        @DisplayName("Checks whether matrix singularity can be determined for edge-size matrices.")
        void doesSingularityWorkForLargeMatrices() {

        }

        @Test
        @DisplayName("Checks whether a matrix is symmetric when off-diagonal elements are nearly the same.")
        void doesSymmetryWorkForCloseValues() {

        }

        @Test
        @DisplayName("Checks whether a matrix is skew symmetric when diagonal elements are nearly 0.")
        void doesSkewSymmetryWorkForNearlyZeroDiagonal() {

        }
    }

    @Nested
    class HelperMethodsTest {

        @Test
        @DisplayName("Test for deep grid copying.")
        void doesDeepCopyWorkAsIntended() {

        }

        @Test
        @DisplayName("Test for bound checking.")
        void doesInBoundsWork() {

        }

        @Test
        @DisplayName("Test for tolerance-based equality.")
        void doesToleranceBasedEqualityWork() {

        }

        @Test
        @DisplayName("Test for row reversal.")
        void doesRowReversalWork() {

        }

        @Test
        @DisplayName("Test for in-place fill.")
        void doesInPlaceFillWorkForLargeMatrix() {

        }

        @Test
        @DisplayName("Test for dot product of vectors.")
        void doesDotProductBehaveCorrectlyForMismatchedVectors() {

        }

        @Test
        @DisplayName("Test for checking jaggedness.")
        void doesIsJaggedThrowForUniformMatrix() {

        }

        @Test
        @DisplayName("Test equalsMatrix method.")
        void doesEqualsMatrixWorkForEdgeSizeMatrices() {

        }

    }

    @Nested
    @DisplayName("Overridden Methods of Object class.")
    class ObjectMethodsTest {

        @Test
        void doNearlyEqualMatricesReturnTrue() {

        }

        @Test
        void doNearlyEqualMatricesHaveSimilarHashcodes() {

        }

        @Test
        void doesToStringHandleEntriesWithDiverseSizes() {

        }

    }
}
