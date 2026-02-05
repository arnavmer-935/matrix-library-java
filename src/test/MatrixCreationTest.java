package test;

import main.Matrix;
import main.MatrixException;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class MatrixCreationTest {

    private static final double TOLERANCE = 1e-6;

    @Test
    @DisplayName("Row/column constructor creates zero-filled matrix with correct dimensions")
    void constructorWithRowsAndColsCreatesZeroMatrix() {
        Matrix A = new Matrix(2,3);

        assertEquals(2, A.getRows());
        assertEquals(3, A.getColumns());

        for (int i = 0; i < A.getRows(); i++) {
            for (int j = 0; j < A.getColumns(); j++) {
                assertEquals(0.0, A.getEntry(i,j), TOLERANCE);
            }
        }
    }

    @Test
    @DisplayName("Square constructor creates matrix with equal row and column count")
    void squareConstructorCreatesCorrectDimensions() {

        Matrix squareMatrix = new Matrix(3);
        assertEquals(squareMatrix.getColumns(), squareMatrix.getRows());
    }

    @Test
    @DisplayName("Constructor rejects negative or zero dimensions")
    void constructorRejectsInvalidDimensions() {

        assertThrows(MatrixException.class, () -> new Matrix(-2,-3));
        assertThrows(MatrixException.class, () -> new Matrix(-5));
        assertThrows(MatrixException.class, () -> new Matrix(0));

    }

    @Test
    @DisplayName("Constructor from 2D array performs deep copy")
    void constructorFrom2DArrayPerformsDeepCopy() {

        double[][] matrixArray = {{1,2,3}, {4,5,6}, {7,8,9}};
        Matrix matrixFrom2dArray = new Matrix(matrixArray);

        matrixArray[0][0] = 67;
        assertEquals(1.0, matrixFrom2dArray.getEntry(0,0), TOLERANCE);
    }

    @Test
    @DisplayName("Constructor rejects null arrays")
    void constructorRejectsNullArrays() {
        double[][] nullGrid = null;
        assertThrows(IllegalArgumentException.class, () -> new Matrix(nullGrid));
    }

    @Test
    @DisplayName("Constructor rejects jagged 2D arrays")
    void constructorRejectsJaggedArrays() {

        double[][] jaggedGrid = {{1,2},{3,4},{5,6,7},{4,5}};
        assertThrows(MatrixException.class, () -> new Matrix(jaggedGrid));
    }

    @Test
    @DisplayName("Constant factory fills matrix with the given value")
    void constantFactoryFillsMatrixWithValue() {

        Matrix constantMatrix = Matrix.constant(4,4,1.2345);
        for (int i = 0; i < constantMatrix.getRows(); i++) {
            for (int j = 0; j < constantMatrix.getColumns(); j++) {
                assertEquals(1.2345, constantMatrix.getEntry(i,j), TOLERANCE);
            }
        }

        Matrix largeConstantMatrix = Matrix.constant(50, 50, -273.6589);
        for (int i = 0; i < constantMatrix.getRows(); i++) {
            for (int j = 0; j < constantMatrix.getColumns(); j++) {
                assertEquals(-273.6589, largeConstantMatrix.getEntry(i,j), TOLERANCE);
            }
        }
    }

    @Test
    @DisplayName("Scalar factory creates 1x1 matrix with the given value")
    void scalarFactoryCreatesUnitMatrix() {

        Matrix scalarMatrix = Matrix.createScalarMatrix(3, 5);
        for (int i = 0; i < scalarMatrix.getRows(); i++) {
            for (int j = 0; j < scalarMatrix.getColumns(); j++) {
                if (i == j) {
                    assertEquals(5, scalarMatrix.getEntry(i,j), TOLERANCE);

                } else {
                    assertEquals(0, scalarMatrix.getEntry(i,j), TOLERANCE);

                }
            }
        }

        Matrix identityMatrix = Matrix.createIdentityMatrix(6);
        for (int i = 0; i < identityMatrix.getRows(); i++) {
            for (int j = 0; j < identityMatrix.getColumns(); j++) {
                if (i == j) {
                    assertEquals(1, identityMatrix.getEntry(i,j), TOLERANCE);

                } else {
                    assertEquals(0, identityMatrix.getEntry(i,j), TOLERANCE);

                }
            }
        }

    }

    @Test
    @DisplayName("ofRows method creates matrix correctly.")
    void ofRowsFactoryMethodCreatesCorrectMatrix() {
        double[] row1 = {1.5, -2, 7.25};
        double[] row2 = {0, 3.14, -999.1};

        double[][] rows = {row1, row2}; //To be used for iterating and testing

        Matrix madeOfRows = Matrix.ofRows(row1, row2);

        assertEquals(2, madeOfRows.getRows());
        assertEquals(3, madeOfRows.getColumns());

        for (int i = 0; i < madeOfRows.getRows(); i++) {
            for (int j = 0; j < madeOfRows.getColumns(); j++) {
                assertEquals(rows[i][j], madeOfRows.getEntry(i,j), TOLERANCE);
            }
        }

        row1[0] = 9999;
        assertEquals(1.5, madeOfRows.getEntry(0,0), TOLERANCE);

    }

    @Test
    @DisplayName("ofColumns method creates matrix correctly.")
    void ofColumnsFactoryMethodCreatesCorrectMatrix() {
        double[] col1 = {1.5, 0};
        double[] col2 = {-2, 3.14};
        double[] col3 = {7.25, -999.1};

        double[] r1 = {1.5, -2, 7.25};
        double[] r2 = {0, 3.14, -999.1};

        double[][] cols = {r1, r2}; //To be used for iterating and testing - represents what the matrix should look like

        Matrix madeOfCols = Matrix.ofColumns(col1, col2, col3);

        assertEquals(2, madeOfCols.getRows());
        assertEquals(3, madeOfCols.getColumns());

        for (int i = 0; i < madeOfCols.getRows(); i++) {
            for (int j = 0; j < madeOfCols.getColumns(); j++) {
                assertEquals(cols[i][j], madeOfCols.getEntry(i,j), TOLERANCE);
            }
        }

        col1[1] = 7253;
        assertEquals(0, madeOfCols.getEntry(1,0), TOLERANCE);
    }
}
