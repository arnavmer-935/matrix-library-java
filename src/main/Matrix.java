package main;
import java.util.*;

/*
* //TODO
* 1. Add test cases for query methods
* 2. Documentation and polish readme
* 3. Tag release
*/

public final class Matrix {

    // ==== INSTANCE VARIABLES AND DECLARATIONS ====
    private record Pair(int x, int y) {
        @Override
        public String toString() {
            return String.format("(%d, %d)", x, y);
        }
    }

    private static class Pivot {

        private static final double PIVOT_TOLERANCE = 1e-12;
        double value;
        int rowIndex;
        int colIndex;
        boolean swapsNeeded;

        Pivot() { //when valid pivot cannot be found
            value = 0.0;
            rowIndex = -1;
            colIndex = -1;
            swapsNeeded = false;
        }

        Pivot(double v, int r, int c, boolean s) {
            value = v;
            rowIndex = r;
            colIndex = c;
            swapsNeeded = s;
        }


        int row() { return rowIndex; }
        boolean swapsNeeded() { return swapsNeeded; }
    }

    private static final double TOLERANCE = 1e-6;
    private final int rows;
    private final int columns;
    private final Pair order;
    private final double[][] entries;


    // ==== MATRIX CREATION METHODS ====
    public Matrix(int nrows, int ncols) {
        if (nrows <= 0 || ncols <= 0) {
            throw MatrixException.illegalDimensions();
        }
        
        this.rows = nrows;
        this.columns = ncols;
        this.order = new Pair(rows, columns);
        this.entries = new double[rows][columns];
    }

    public Matrix(double[][] grid) throws MatrixException {

        validateGrid(grid);
        this.entries = deepGridCopy(grid);
        this.rows = grid.length;
        this.columns = grid[0].length;
        this.order = new Pair(rows, columns);
    }

    public Matrix(List<List<Double>> grid) {

        if (grid == null) {
            throw new IllegalArgumentException("Matrix grid must be non-null.");
        }

        if (grid.isEmpty() || grid.getFirst().isEmpty()) {
            throw MatrixException.illegalDimensions();
        }

        if (isJaggedGrid(grid)) {
            throw MatrixException.jaggedMatrix(getMismatchedRowIndex(grid));
        }

        this.rows = grid.size();
        this.columns = grid.getFirst().size();
        this.order = new Pair(rows, columns);
        this.entries = new double[rows][columns];

        for (int i = 0; i < rows; i++) {
            List<Double> ithRow = grid.get(i);
            if (ithRow == null) {
                throw new IllegalArgumentException("Matrix grid must be non-null.");
            }

            for (int j = 0; j < columns; j++) {
                if (!Double.isFinite(ithRow.get(j))) {
                    throw new IllegalArgumentException("Infinite/undefined value at coordinates " + new Pair(i,j));
                }
                
                this.entries[i][j] = ithRow.get(j);
            }
        }
    }

    //Constructor for square matrix
    public Matrix(int nrows) {
        if (nrows <= 0) {
            throw MatrixException.illegalDimensions();
        }
        this.rows = nrows;
        this.columns = nrows;
        this.order = new Pair(rows, rows);
        this.entries = new double[rows][columns];
    }

    //deep copy constructor
    public Matrix(Matrix A) {
        this(getValidMatrix(A));
    }

    public static Matrix ofRows(double[]... rows) {
        validateGrid(rows);
        double[][] result = new double[rows.length][rows[0].length];

        for (int i = 0; i < rows.length; i++) {
            for (int j = 0; j < rows[0].length; j++) {
                result[i][j] = rows[i][j];
            }
        }

        return new Matrix(result);
    }

    public static Matrix ofColumns(double[]... columns) {
        return Matrix.ofRows(columns).transpose();
    }

    public static Matrix constant(int nrows, int ncols, double k) {
        if (nrows <= 0 || ncols <= 0) {
            throw MatrixException.illegalDimensions();
        }

        if (!Double.isFinite(k)) {
            throw MatrixException.infiniteValue();
        }

        if (almostEqual(k,0)) {
            return new Matrix(nrows, ncols);
        }

        double[][] res = new double[nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            Arrays.fill(res[i], k);
        }

        return new Matrix(res);
    }

    public void zeroMatrixInPlace() {
        fillInPlace(0.0);
    }

    public static Matrix createScalarMatrix(int nrows, double k) {
        if (!Double.isFinite(k)) {
            throw MatrixException.infiniteValue();
        }

        if (nrows <= 0) {
            throw MatrixException.illegalDimensions();
        }

        if (almostEqual(k,0)) {
            return new Matrix(nrows);
        }

        double[][] result = new double[nrows][nrows];
        for (int i = 0; i < nrows; i++) {
            result[i][i] = k;
        }

        return new Matrix(result);
    }

    public static Matrix createIdentityMatrix(int nrows) {
        return createScalarMatrix(nrows, 1.0);
    }


    // ==== ACCESSOR AND MUTATOR METHODS ====
    public int getRows() {
        return rows;
    }

    public int getColumns() {
        return columns;
    }

    public Pair getOrder() { return this.order; }

    private double[][] getEntries() { return this.entries; }

    public double getEntry(int rowIndex, int colIndex) {
        return getValue(rowIndex, colIndex);
    }

    public void setEntry(double value, int r, int c) {
        if (!isInBounds(r,c)) {
            throw new IndexOutOfBoundsException(String.format("Cell coordinates (%d, %d) out of bounds for Matrix of order %s", r, c, order));
        }

        if (!Double.isFinite(value)) {
            throw new IllegalArgumentException("Given value cannot be used in matrix since it is Not a Number (Infinite/undefined value)");
        }

        this.entries[r][c] = value;
    }

    public void setRow(double[] row, int rowIndex) {

        if (row == null) {
            throw new IllegalArgumentException("Input row must be non-null.");
        }

        if (row.length != columns) {
            throw MatrixException.rowLengthMismatch();
        }

        if (!rowInRange(rowIndex)) {
            throw new IndexOutOfBoundsException(String.format("Row index %d is out of bounds for Matrix of order %s", rowIndex, this.order));
        }

        for (int i = 0; i < row.length; i++) {
            if (!Double.isFinite(row[i])) {
                throw new IllegalArgumentException("Infinite/undefined value found at row index " + i);
            }
        }

        this.entries[rowIndex] = row.clone();
    }

    public void setColumn(double[] col, int colIndex) {
        if (col == null) {
            throw new IllegalArgumentException("Input column must be non-null.");
        }

        if (col.length != rows) {
            throw MatrixException.columnLengthMismatch();

        }

        if (!colInRange(colIndex)) {
            throw new IndexOutOfBoundsException(String.format("Column index %d is out of bounds for Matrix of order %s", colIndex, this.order));
        }

        for (int i = 0; i < col.length; i++) {
            if (!Double.isFinite(col[i])) {
                throw new IllegalArgumentException("Infinite/undefined value found at row index " + i);
            }
        }

        for (int i = 0; i < rows; i++) {
            this.entries[i][colIndex] = col[i];
        }
    }

    // ==== IN-PLACE OPERATIONS ====
    public void multiplyByScalarInPlace(double k) {
        if (!Double.isFinite(k)) {
            throw MatrixException.infiniteValue();
        }

        if (almostEqual(k, 1.0)) {
            return;
        }

        if (almostEqual(k,0.0)) {
            zeroMatrixInPlace();
            return;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                entries[i][j] *= k;
            }
        }
    }

    public void addInPlace(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operands must be non-null.");
        }

        if (!this.order.equals(other.order)) {
            throw MatrixException.orderMismatch();
        }

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {
                this.entries[r][c] += other.entries[r][c];
            }
        }
    }

    public void subtractInPlace(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operands must be non-null.");
        }

        if (!this.order.equals(other.order)) {
            throw MatrixException.orderMismatch();
        }

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {
                this.entries[r][c] -= other.entries[r][c];
            }
        }
    }

    public void transposeInPlace() {
        if (!isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        for (int i = 0; i < rows; i++) {
            for (int j = i + 1; j < columns; j++) {
                double temp = entries[i][j];
                entries[i][j] = entries[j][i];
                entries[j][i] = temp;
            }
        }
    }

    // ==== OUT OF PLACE OPERATIONS ====
    public Matrix multiplyByScalar(double k) {
        if (!Double.isFinite(k)) {
            throw MatrixException.infiniteValue();
        }

        if (almostEqual(k,0.0)) {
            return new Matrix(rows, columns);
        }

        if (almostEqual(k,1.0)) {
            return new Matrix(entries);
        }

        double[][] result = new double[this.rows][this.columns];
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = (this.entries[i][j] * k);
            }
        }
        return new Matrix(result);
    }

    public Matrix add(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operand must be non-null.");
        }

        if (!this.order.equals(other.order)) {
            throw MatrixException.orderMismatch();
        }

        double[][] result = new double[other.rows][other.columns];
        for (int r = 0; r < result.length; r++) {
            for (int c = 0; c < result[0].length; c++) {
                result[r][c] = (this.entries[r][c] + other.entries[r][c]);
            }
        }
        return new Matrix(result);
    }

    public Matrix subtract(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operand must be non-null.");
        }
        return add(other.multiplyByScalar(-1));
    }

    public Matrix transpose() { //used for transposing rectangular matrices
        double[][] transposedGrid = new double[this.columns][this.rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                transposedGrid[j][i] = this.entries[i][j];
            }
        }
        return new Matrix(transposedGrid);
    }

    public Matrix multiply(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operand must be non-null.");
        }

        if (this.columns != other.rows) {
            throw MatrixException.dimensionMismatch();
        }

        if (this.isZeroMatrix() || other.isZeroMatrix()) {
            return new Matrix(this.rows, other.columns);
        }

        double[][] product = new double[this.rows][other.columns];

        for (int j = 0; j < other.columns; j++) {
            double[] jthColumn = other.getColumn(j);
            for (int i = 0; i < this.rows; i++) {
                product[i][j] = dotProduct(this.entries[i], jthColumn);
            }
        }

        return new Matrix(product);
    }

    public double[][] toArray() {
        return deepGridCopy(this.entries);
    }

    // ==== NUMERICAL METHODS ====

    public Matrix inverse() {

        if (!isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        double[][] left = deepGridCopy(this.entries);
        double[][] right = createIdentityMatrix(rows).toArray();

        for (int i = 0; i < left[0].length; i++) {
            Pivot pivot = findValidPivot(left, i);

            if (pivot.row() == -1) {
                throw MatrixException.matrixSingularity();
            }

            if (pivot.swapsNeeded()) {
                swapGridRow(left, i, pivot.row());
                swapGridRow(right, i, pivot.row());
            }

            double scaleFactor = left[i][i];
            for (int c = 0; c < left[0].length; c++) { //normalizing pivot row
                left[i][c] /= scaleFactor;
                right[i][c] /= scaleFactor;
            }

            for (int row = 0; row < left.length; row++) {
                if (row == i) continue;
                double factor = left[row][i]; //since pivot value becomes 1 after normalization;

                for (int j = 0; j < left[row].length; j++) {
                    left[row][j] -= (factor * left[i][j]);
                    right[row][j] -= (factor * right[i][j]);
                }
            }
        }

        return new Matrix(right);

    }

    public double determinant() {
        if (!isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        if (rows == 1 && columns == 1) { //single element matrix
            return entries[0][0];
        }

        if (rows == 2 && columns == 2) { //2x2 square matrix
            return (entries[0][0] * entries[1][1]) - (entries[0][1] * entries[1][0]);
        }

        if (isUpperTriangular()) {
            return diagonalProduct(this.entries);
        }

        double[][] grid = deepGridCopy(this.entries); //copying to preserve original
        int ROWS = grid.length;
        int COLS = grid[0].length;
        int swapCount = 0;

        for (int c = 0; c < COLS; c++) {

            Pivot pivot = findValidPivot(grid, c);

            if (pivot.row() == -1) {
                return 0.0;
            }

            if (pivot.swapsNeeded()) {
                swapGridRow(grid, c, pivot.row());
                swapCount++;
            }

            for (int rowCursor = c + 1; rowCursor < ROWS; rowCursor++) {
                double factor = grid[rowCursor][c] / grid[c][c];
                for (int i = 0; i < grid[rowCursor].length; i++) {

                    double elementInPivotRow = grid[c][i];
                    grid[rowCursor][i] -= (factor * elementInPivotRow);
                }
            }
        }
        double sign = swapCount % 2 == 0 ? 1.0: -1.0;
        return sign * diagonalProduct(grid);

    }

    // ==== QUERY METHODS ====
    public boolean isSquareMatrix() {
        return this.rows == this.columns;
    }

    public boolean isIdentityMatrix() {
        if (!isSquareMatrix()) {
            return false;
        }

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {

                if (r != c && !almostEqual(entries[r][c], 0)) return false;
                if (r == c && !almostEqual(entries[r][c], 1.0)) return false;

            }
        }
        return true;
    }

    public boolean isUpperTriangular() {
        if (!isSquareMatrix()) {
            return false;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < i; j++) {
                if (!almostEqual(entries[i][j], 0.0)) return false;
            }
        }

        return true;
    }

    public boolean isLowerTriangular() {
        if (!isSquareMatrix()) {
            return false;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = i+1; j < columns; j++) {
                if (!almostEqual(entries[i][j], 0.0)) return false;
            }
        }

        return true;
    }

    public boolean isZeroMatrix() {
        return areAllEqual(0.0);
    }

    public boolean isConstantMatrix() {
        return areAllEqual(entries[0][0]);
    }

    public boolean isSingular() {
        return almostEqual(determinant(), 0.0);
    }

    public boolean isSymmetric() {

        if (!isSquareMatrix()) {
            return false;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = i+1; j < columns; j++) {
                if (!almostEqual(entries[i][j], entries[j][i])) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isSkewSymmetric() {

        if (!isSquareMatrix()) {
            return false;
        }

        if (!onlyZeroesInDiagonal()) {
            return false;
        }

        for (int i = 0; i < this.rows; i++) {
            for (int j = i+1; j < this.columns; j++) {
                if (!almostEqual(entries[i][j], -entries[j][i])) {
                    return false;
                }
            }
        }
        return true;
    }

    // ==== HELPER METHODS ====
    private boolean isInBounds(int r, int c) {
        return rowInRange(r) && colInRange(c);
    }

    private double getValue(int r, int c) {
        if (!isInBounds(r, c)) {
            throw new IndexOutOfBoundsException(String.format("Cell coordinates (%d, %d) out of bounds for Matrix of order %s", r, c, order));
        }
        return this.entries[r][c];
    }

    private static boolean almostEqual(double a, double b) {
        return Math.abs(a-b) <= TOLERANCE;
    }

    private static void reverseRow(double[] row) {
        int n = row.length;
        for (int i = 0; i < n / 2; i++) {
            double temp = row[i];
            row[i] = row[n - i - 1];
            row[n - i - 1] = temp;
        }
    }

    private void fillInPlace(double k) {
        for (int i = 0; i < this.rows; i++) {
            Arrays.fill(entries[i], k);
        }
    }

    private boolean areAllEqual(double k) {
        for (double[] row : this.entries) {
            for (double n : row) {
                if (!almostEqual(n, k)) {
                    return false;
                }
            }
        }
        return true;
    }

    private static double dotProduct(double[] arr1, double[] arr2) {
        int len1 = arr1.length;
        int len2 = arr2.length;
        double res = 0;
        for (int i = 0; i < len1; i++) {
            res += (arr1[i] * arr2[i]);
        }
        return res;
    }

    private double[] getColumn(int n) {

        double[] col = new double[this.rows];
        for (int i = 0; i < this.rows; i++) {
            col[i] = this.entries[i][n];
        }
        return col;
    }

    private static boolean isJaggedGrid(List<List<Double>> grid) {
        return getMismatchedRowIndex(grid) != -1;
    }

    private static boolean isJaggedGrid(double[][] grid) {
        return getMismatchedRowIndex(grid) != -1;
    }


    private static int getMismatchedRowIndex(double[][] grid) {
        int refSize = grid[0].length;
        for (int i = 1; i < grid.length; i++) {
            if (grid[i].length != refSize) {
                return i;
            }
        }
        return -1; //No such row found
    }

    private static int getMismatchedRowIndex(List<List<Double>> grid) {
        int refSize = grid.getFirst().size();
        for (int i = 1; i < grid.size(); i++) {
            if (grid.get(i).size() != refSize) {
                return i;
            }
        }
        return -1; //No such row whose length is different from reference size
    }

    private boolean onlyZeroesInDiagonal() {
        for (int i = 0; i < this.rows; i++) {
            if (!almostEqual(entries[i][i], 0.0)) {
                return false;
            }
        }
        return true;
    }

    private boolean equalsMatrix(Matrix other) {
        if (this == other) return true;
        if (other == null || !this.order.equals(other.order)) return false;

        double[][] tempGrid = other.getEntries();
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                if (Double.compare(this.entries[i][j], tempGrid[i][j]) != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    private boolean rowInRange(int rowIndex) {
        return 0 <= rowIndex && rowIndex < rows;
    }

    private boolean colInRange(int colIndex) {
        return 0 <= colIndex && colIndex < columns;
    }

    private double[][] deepGridCopy(double[][] grid) {

        int nrows = grid.length;
        int ncols = grid[0].length;
        double[][] result = new double[nrows][ncols];

        for (int i = 0; i < nrows; i++) {
            System.arraycopy(grid[i], 0, result[i], 0, ncols);
        }
        return result;
    }

    private Pivot findValidPivot(double[][] grid, int colIdx) {

        if (Math.abs(grid[colIdx][colIdx]) > Pivot.PIVOT_TOLERANCE) {
            return new Pivot(grid[colIdx][colIdx], colIdx, colIdx, false);
        }

        for (int rowIdx = colIdx + 1; rowIdx < grid.length; rowIdx++) {
            double possibleValue = grid[rowIdx][colIdx];
            if (!almostEqual(possibleValue, 0.0)) {
                return new Pivot(possibleValue, rowIdx, colIdx, true);
            }
        }
        return new Pivot();
    }

    private void swapGridRow(double[][] grid, int c, int row) {
        double[] temp = grid[c];
        grid[c] = grid[row];
        grid[row] = temp;
    }

    private static double diagonalProduct(double[][] grid) {
        double res = 1;
        for (int i = 0; i < grid.length; i++) {
            res *= grid[i][i];
        }
        return res;
    }

    private String formattedValue(double v) {
        return Math.abs(v) <= TOLERANCE ? "0.0000" : String.format("%.4f", v);
    }

    private static boolean containsNullRows(double[][] grid) {
        for (double[] row : grid) {
            if (row == null) {
                return true;
            }
        }
        return false;
    }

    private static void validateGrid(double[][] grid) {
        if (grid == null || containsNullRows(grid)) {
            throw new IllegalArgumentException("Matrix grid must be non-null.");
        }

        if (grid.length == 0 || grid[0].length == 0) {
            throw MatrixException.illegalDimensions();
        }

        if (isJaggedGrid(grid)) {
            throw MatrixException.jaggedMatrix(getMismatchedRowIndex(grid));
        }

        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[0].length; j++) {
                if (!Double.isFinite(grid[i][j])) {
                    throw new IllegalArgumentException("Infinite value at coordinates " + new Pair(i,j));
                }
            }
        }
    }

    private static double[][] getValidMatrix(Matrix A) {
        if (A == null) {
            throw new IllegalArgumentException("Matrix input must be non-null.");
        }

        return A.toArray();
    }

    // ==== OBJECT METHODS ====

    @Override
    public String toString() {

        String[][] formattedValues = new String[rows][columns];
        int longestLength = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                String val = formattedValue(entries[i][j]);
                longestLength = Math.max(longestLength, val.length());

                formattedValues[i][j] = val;
            }
        }

        StringBuilder res = new StringBuilder();
        for (int r = 0; r < rows; r++) {
            res.append("[ ");
            for (int c = 0; c < columns; c++) {
                res.append(String.format("%" + longestLength + "s", formattedValues[r][c]));
                if (c < columns - 1) res.append("  ");
            }
            res.append(" ]\n");
        }
        return res.toString();
    }

    @Override
    public int hashCode() {
        return 31 * order.hashCode() + Arrays.deepHashCode(this.entries);
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof Matrix)) return false;

        Matrix o = (Matrix)other;
        return this.equalsMatrix(o);
    }

}

