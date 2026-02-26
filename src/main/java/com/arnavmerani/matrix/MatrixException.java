package main.java.com.arnavmerani.matrix;
/**
 * An unchecked exception indicating that a matrix operation
 * could not be completed due to invalid structure, dimensions,
 * or numerical conditions.
 * <p>
 * Contains factory methods for throwing exceptions, with the aim
 * to promote readability and reusability.
 *
 * <p>This exception is thrown in cases such as:
 * <ul>
 *   <li>Uneven row lengths during initialization</li>
 *   <li>Dimension mismatches in arithmetic operations</li>
 *   <li>Operations requiring a square matrix applied to a non-square matrix</li>
 *   <li>Attempts to invert a singular matrix</li>
 * </ul>
 *
 * <p>{@code MatrixException} extends {@code RuntimeException}
 * and is intended to signal programming errors or invalid usage
 * of the {@code Matrix} API.
 */
public class MatrixException extends RuntimeException {
    private MatrixException(String message) {
        super(message);
    }

    public static MatrixException illegalDimensions() {
        return new MatrixException("Matrix dimensions must be positive.");
    }

    public static MatrixException jaggedMatrix(int mismatchedIndex) {
        return new MatrixException(String.format("Length of row at index %d does not match the length of row at index 0",
                mismatchedIndex));
    }

    public static MatrixException rowLengthMismatch() {
        return new MatrixException("Row length does not match the number of columns in the matrix.");
    }

    public static MatrixException columnLengthMismatch() {
        return new MatrixException("Column length does not match the number of rows in the matrix.");
    }

    public static MatrixException orderMismatch() {
        return new MatrixException("The orders of the given operands must match.");
    }

    public static MatrixException requireSquareMatrix() {
        return new MatrixException("This operation requires the Matrix to have an equal number of rows and columns.");
    }

    public static MatrixException dimensionMismatch() {
        return new MatrixException("Number of columns in first matrix must match the number of rows in second matrix.");
    }

    public static MatrixException matrixSingularity() {
        return new MatrixException("This operation is not defined for a singular matrix.");
    }

    public static MatrixException infiniteValue() {
        return new MatrixException("Scalar must be a finite value.");
    }
}