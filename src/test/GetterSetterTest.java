package test;

import main.Matrix;
import main.MatrixException;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static java.lang.Double.NaN;
import static org.junit.jupiter.api.Assertions.*;

public class GetterSetterTest {

    @Test
    @DisplayName("setEntry method throws a MatrixException when trying to update a matrix entry.")
    void setEntryThrowsExceptionForOutOfBounds() {
        assertThrows(IndexOutOfBoundsException.class, () -> new Matrix(4,5).setEntry(5.0, 4,4));
        assertThrows(IndexOutOfBoundsException.class, () -> new Matrix(4,5).setEntry(6.0, 3,5));

        assertThrows(IndexOutOfBoundsException.class, () -> new Matrix(1,1).setEntry(5.0, 1,1));
        assertThrows(IndexOutOfBoundsException.class, () -> new Matrix(50,50).setEntry(4.99, 49,50));

        Matrix invalidCoordinates = new Matrix(3,3);
        assertThrows(IndexOutOfBoundsException.class, () -> invalidCoordinates.setEntry(2.0, -1,0));
        assertThrows(IndexOutOfBoundsException.class, () -> invalidCoordinates.setEntry(4.0, 1,-2));
        assertThrows(IndexOutOfBoundsException.class, () -> invalidCoordinates.setEntry(4.0, -1,-1));
    }

    @Test
    @DisplayName("setRow works")
    void setRowUpdatesMatrixCorrectly() {
        double[] newRow = {1.0, 3.5, -19.5, 0.0};
        double[] invalidRow = {1.0, 0.0, 9.2, NaN};
        double[] exceptionRow = {1.0, 3.0, -4.1};

        Matrix setRowTest = new Matrix(3,4);

        setRowTest.setRow(newRow, 0);
        newRow[0] = 6.7;
        assertEquals(1.0, setRowTest.getEntry(0,0)); //tests copying

        assertThrows(IllegalArgumentException.class, () -> new Matrix(4,4).setRow(invalidRow, 0));
        assertThrows(MatrixException.class, () -> new Matrix(3,2).setRow(exceptionRow, 0));
        assertThrows(IndexOutOfBoundsException.class, () -> new Matrix(3,3).setRow(exceptionRow, 3));


    }

    @Test
    @DisplayName("setCol works")
    void setColUpdatesMatrixCorrectly() {
        double[] newCol = {1.0, 3.5, -19.5, 0.0};
        double[] invalidCol = {1.0, 0.0, 9.2, NaN};
        double[] exceptionCol = {1.0, 3.0, -4.1};

        Matrix setColTest = new Matrix(4,3);

        setColTest.setColumn(newCol, 0);
        newCol[0] = 6.7;
        assertEquals(1.0, setColTest.getEntry(0,0)); //tests copying

        assertThrows(IllegalArgumentException.class, () -> new Matrix(4,4).setColumn(invalidCol, 0));
        assertThrows(MatrixException.class, () -> new Matrix(2,3).setColumn(exceptionCol, 0));
        assertThrows(IndexOutOfBoundsException.class, () -> new Matrix(3,3).setColumn(exceptionCol, 3));

    }

    @Test
    @DisplayName("getEntry works with bounds checks")
    void getEntryReturnsCorrectly() {
        double[] r1 = {1.0, 3.5, -19.5, 0.0};
        double[] r2 = {1.0, 0.0, 9.2, 6.78};
        double[] r3 = {1.0, 3.0, -4.1, 9.385};

        Matrix testGetEntry = Matrix.ofRows(r1, r2, r3);

        assertEquals(r1[0], testGetEntry.getEntry(0,0));
        assertEquals(r2[2], testGetEntry.getEntry(1,2));
        assertThrows(IndexOutOfBoundsException.class, () -> testGetEntry.getEntry(3,4));
        assertThrows(IndexOutOfBoundsException.class, () -> testGetEntry.getEntry(3,2));
        assertThrows(IndexOutOfBoundsException.class, () -> testGetEntry.getEntry(1,4));
        assertThrows(IndexOutOfBoundsException.class, () -> testGetEntry.getEntry(-1,-2));
        assertThrows(IndexOutOfBoundsException.class, () -> testGetEntry.getEntry(-2,0));
    }
}
