package test;

import main.Matrix;
import main.MatrixException;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static java.lang.Double.NaN;
import static org.junit.jupiter.api.Assertions.*;

public class GetterSetterTest {

    @Nested
    class GetEntryTests {

        @Test
        @DisplayName("getEntry returns correct value for valid indices")
        void getEntry_validIndices() {
            Matrix M = Matrix.ofRows(
                    new double[]{1.0, 3.5, -19.5, 0.0},
                    new double[]{1.0, 0.0, 9.2, 6.78},
                    new double[]{1.0, 3.0, -4.1, 9.385}
            );

            assertEquals(1.0, M.getEntry(0,0));
            assertEquals(9.2, M.getEntry(1,2));
        }

        @Test
        @DisplayName("getEntry throws for invalid indices")
        void getEntry_invalidIndices() {
            Matrix M = new Matrix(3,4);

            assertThrows(IndexOutOfBoundsException.class, () -> M.getEntry(3,4));
            assertThrows(IndexOutOfBoundsException.class, () -> M.getEntry(3,2));
            assertThrows(IndexOutOfBoundsException.class, () -> M.getEntry(1,4));
            assertThrows(IndexOutOfBoundsException.class, () -> M.getEntry(-1,-2));
            assertThrows(IndexOutOfBoundsException.class, () -> M.getEntry(-2,0));
        }
    }

    @Nested
    class SetEntryTests {

        @Test
        @DisplayName("setEntry updates correct position")
        void setEntry_valid() {
            Matrix M = new Matrix(2,2);
            M.setEntry(5.0, 1, 1);
            assertEquals(5.0, M.getEntry(1,1));
        }

        @Test
        @DisplayName("setEntry throws for out-of-bounds indices")
        void setEntry_invalidIndices() {
            Matrix M = new Matrix(4,5);

            assertThrows(IndexOutOfBoundsException.class, () -> M.setEntry(5.0, 4,4));
            assertThrows(IndexOutOfBoundsException.class, () -> M.setEntry(6.0, 3,5));

            Matrix single = new Matrix(1,1);
            assertThrows(IndexOutOfBoundsException.class, () -> single.setEntry(5.0, 1,1));

            Matrix invalid = new Matrix(3,3);
            assertThrows(IndexOutOfBoundsException.class, () -> invalid.setEntry(2.0, -1,0));
            assertThrows(IndexOutOfBoundsException.class, () -> invalid.setEntry(4.0, 1,-2));
            assertThrows(IndexOutOfBoundsException.class, () -> invalid.setEntry(4.0, -1,-1));
        }
    }

    @Nested
    class SetRowTests {

        @Test
        @DisplayName("setRow updates matrix and performs defensive copy")
        void setRow_validAndCopies() {
            Matrix M = new Matrix(3,4);
            double[] row = {1.0, 3.5, -19.5, 0.0};

            M.setRow(row, 0);
            row[0] = 999.0;

            assertEquals(1.0, M.getEntry(0,0));
        }

        @Test
        @DisplayName("setRow throws for invalid input")
        void setRow_invalidInput() {
            double[] invalidRow = {1.0, 0.0, 9.2, NaN};
            double[] wrongLength = {1.0, 3.0, -4.1};

            assertThrows(IllegalArgumentException.class,
                    () -> new Matrix(4,4).setRow(invalidRow, 0));

            assertThrows(MatrixException.class,
                    () -> new Matrix(3,2).setRow(wrongLength, 0));

            assertThrows(IndexOutOfBoundsException.class,
                    () -> new Matrix(3,3).setRow(wrongLength, 3));

            assertThrows(IllegalArgumentException.class,
                    () -> new Matrix(2,3).setRow(null, 0));
        }

        @Test
        @DisplayName("setRow is atomic on failure")
        void setRow_atomicity() {
            Matrix M = new Matrix(3,4);
            Matrix copy = new Matrix(M.toArray());
            double[] badRow = {1, -2, 3, NaN};

            assertThrows(IllegalArgumentException.class,
                    () -> M.setRow(badRow, 0));

            assertEquals(copy, M);
        }
    }

    @Nested
    class SetColumnTests {

        @Test
        @DisplayName("setColumn updates matrix and performs defensive copy")
        void setColumn_validAndCopies() {
            Matrix M = new Matrix(4,3);
            double[] col = {1.0, 3.5, -19.5, 0.0};

            M.setColumn(col, 0);
            col[0] = 999.0;

            assertEquals(1.0, M.getEntry(0,0));
        }

        @Test
        @DisplayName("setColumn throws for invalid input")
        void setColumn_invalidInput() {
            double[] invalidCol = {1.0, 0.0, 9.2, NaN};
            double[] wrongLength = {1.0, 3.0, -4.1};

            assertThrows(IllegalArgumentException.class,
                    () -> new Matrix(4,4).setColumn(invalidCol, 0));

            assertThrows(MatrixException.class,
                    () -> new Matrix(2,3).setColumn(wrongLength, 0));

            assertThrows(IndexOutOfBoundsException.class,
                    () -> new Matrix(3,3).setColumn(wrongLength, 3));

            assertThrows(IllegalArgumentException.class,
                    () -> new Matrix(3,3).setColumn(null, 1));
        }

        @Test
        @DisplayName("setColumn is atomic on failure")
        void setColumn_atomicity() {
            Matrix M = new Matrix(3,4);
            Matrix copy = new Matrix(M.toArray());
            double[] badCol = {1, -2, NaN};

            assertThrows(IllegalArgumentException.class,
                    () -> M.setColumn(badCol, 0));

            assertEquals(copy, M);
        }
    }
}