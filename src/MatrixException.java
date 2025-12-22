public class MatrixException extends RuntimeException {
    private MatrixException(String message) {
        super(message);
    }

    public static MatrixException illegalDimensions() {
        return new MatrixException("Matrix dimensions must be positive.");
    }

    public static MatrixException jaggedMatrix(int mismatchedIndex) {
        return new MatrixException(String.format("Length of row at index %d does not match the length of row at index 0"));
    }
}