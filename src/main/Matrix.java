package main;
import java.util.*;

public class Matrix {

    // ==== INSTANCE VARIABLES AND DECLARATIONS ====
    private record Pair(int x, int y) {
        @Override
        public String toString() {
            return String.format("(%d, %d)", x, y);
        }
    }

    private static class Pivot {
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
    private double[][] entries;


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

        this.entries = deepGridCopy(grid);
        this.rows = grid.length;
        this.columns = grid[0].length;
        this.order = new Pair(rows, columns);
    }

    public Matrix(List<List<Double>> grid) {

        if (grid == null || containsNullRows(grid)) {
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

    public static Matrix ofRows(double[]... rows) {
        if (rows == null || containsNullRows(rows)) {
            throw new IllegalArgumentException("Matrix grid must be non-null.");
        }

        if (rows.length == 0 || rows[0].length == 0) {
            throw MatrixException.illegalDimensions();
        }

        if (isJaggedGrid(rows)) {
            throw MatrixException.jaggedMatrix(getMismatchedRowIndex(rows));
        }
        
        double[][] result = new double[rows.length][rows[0].length];
        for (int i = 0; i < rows.length; i++) {
            for (int j = 0; j < rows[0].length; j++) {
                if (!Double.isFinite(rows[i][j])) {
                    throw MatrixException.infiniteValue();
                }
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

        double[][] res = new double[nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            Arrays.fill(res[i], k);
        }

        return new Matrix(res);
    }

    public static Matrix zeroMatrix(int nrows, int ncols) {
        return constant(nrows, ncols, 0.0);
    }

    public void nullMatrixInPlace() {
        fillInPlace(0.0);
    }

    public static Matrix createScalarMatrix(int nrows, double k) {
        if (!Double.isFinite(k)) {
            throw MatrixException.infiniteValue();
        }

        if (nrows <= 0) {
            throw MatrixException.illegalDimensions();
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

        if (k == 1.0) {
            return;
        }

        if (k == 0.0) {
            this.nullMatrixInPlace();
            return;
        }

        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                entries[i][j] *= k;
            }
        }
    }

    public void addInPlace(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operands must be non-null.");
        }

        if (!this.getOrder().equals(other.getOrder())) {
            throw MatrixException.orderMismatch();
        }

        for (int r = 0; r < this.rows; r++) {
            for (int c = 0; c < this.columns; c++) {
                this.entries[r][c] += other.getValue(r, c);
            }
        }
    }

    public void subtractInPlace(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operands must be non-null.");
        }

        if (!this.getOrder().equals(other.getOrder())) {
            throw MatrixException.orderMismatch();
        }

        for (int r = 0; r < this.rows; r++) {
            for (int c = 0; c < this.columns; c++) {
                this.entries[r][c] -= other.getValue(r, c);
            }
        }
    }

    public void transposeInPlace() {
        if (!this.isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        for (int i = 0; i < rows; i++) {
            for (int j = i+1; j < columns; j++) {
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

        if (k == 0.0) {
            return zeroMatrix(this.rows, this.columns);
        }

        if (k == 1.0) {
            return new Matrix(this.entries);
        }

        double[][] result = new double[this.rows][this.columns];
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = this.getValue(i,j) * k;
            }
        }
        return new Matrix(result);
    }

    public Matrix add(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operand must be non-null.");
        }

        if (!this.getOrder().equals(other.getOrder())) {
            throw MatrixException.orderMismatch();
        }

        double[][] result = new double[other.rows][other.columns];
        for (int r = 0; r < result.length; r++) {
            for (int c = 0; c < result[0].length; c++) {
                result[r][c] = this.getValue(r, c) + other.getValue(r, c);
            }
        }
        return new Matrix(result);
    }

    public Matrix subtract(Matrix other) {
        return add(other.multiplyByScalar(-1));
    }

    public Matrix transpose() { //used for transposing rectangular matrices
        double[][] transposedGrid = new double[this.columns][this.rows];

        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
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

        if (this.isNullMatrix() || other.isNullMatrix()) {
            return zeroMatrix(this.rows, other.columns);
        }

        double[][] product = new double[this.rows][other.columns];

        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                double[] jthColumn = getColumn(j);
                product[i][j] = dotProduct(this.entries[i], jthColumn);
            }
        }

        return new Matrix(product);
    }

    public Matrix inverse() {

        if (isJaggedGrid(this.entries)) {
            throw MatrixException.jaggedMatrix(getMismatchedRowIndex(this.entries));
        }

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

    public double[][] toArray() {
        return deepGridCopy(this.entries);
    }

    // ==== NUMERICAL METHODS ====

    public double determinant() {
        if (isJaggedGrid(this.entries)) {
            throw MatrixException.jaggedMatrix(getMismatchedRowIndex(this.entries));
        }

        if (!this.isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        if (this.rows == 1 && this.columns == 1) { //single element matrix
            return this.entries[0][0];
        }

        if (this.rows == 2 && this.columns == 2) { //2x2 square matrix
            return (this.entries[0][0] * this.entries[1][1]) - (this.entries[0][1] * this.entries[1][0]);
        }

        double[][] grid = deepGridCopy(this.entries); //copying to preserve original
        int ROWS = grid.length;
        int COLS = grid[0].length;
        int swapCount = 0;

        for (int c = 0; c < COLS; c++) {
            //first find a valid pivot
            Pivot pivot = findValidPivot(grid, c);

            if (pivot.row() == -1) { //valid, non-zero pivot not found in column, hence the determiInfinite/undefined valuet is zero
                return 0.0;
            }

            if (pivot.swapsNeeded()) {
                swapGridRow(grid, c, pivot.row());
                swapCount++;
            }

            //if it lies on diagonal, no swaps needed (covered in findValidPivot method)
            // otherwise, swap pivot row with the row at index 'c'.
            //increment swap count

            //then iterate over subsequent elements of column c.
            for (int rowCursor = c + 1; rowCursor < ROWS; rowCursor++) {
                //calculate row reduction factor and apply to entire row
                for (int i = 0; i < grid[rowCursor].length; i++) {

                    double factor = grid[rowCursor][c] / grid[c][c];
                    double elementInPivotRow = grid[c][i];
                    grid[rowCursor][i] -= (factor * elementInPivotRow);
                }
            }
        }

        return Math.pow(-1, swapCount) * diagonalProduct(grid);
    }

    // ==== QUERY METHODS ====
    public boolean isSquareMatrix() {
        return this.rows == this.columns;
    }

    public boolean isIdentityMatrix() {
        if (!this.isSquareMatrix()) {
            return false;
        }

        for (int i = 0; i < this.rows; ++i) {
            if (!almostEqual(entries[i][i], 1.0)) {
                return false;
            }
        }

        for (int r = 0; r < this.rows; r++) {
            for (int c = 0; c < this.columns; c++) {
                if (r != c && !almostEqual(entries[r][c], 0)) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isUpperTriangular() {
        for (Double x : this.getLowerTriangle()) {
            if ((!almostEqual(x, 0.0))) {
                return false;
            }
        }
        return true;
    }

    public boolean isLowerTriangular() {
        for (Double x : this.getUpperTriangle()) {
            if (!almostEqual(x,0.0)) {
                return false;
            }
        }
        return true;
    }

    public boolean isNullMatrix() {
        return this.areAllEqual(0.0);
    }

    public boolean isConstantMatrix() {
        return this.areAllEqual(entries[0][0]);
    }

    public boolean isSingular() {
        return almostEqual(determinant(), 0.0);
    }

    public boolean preceeds(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operand must be non-null.");
        }

        if (!this.getOrder().equals(other.getOrder())) {
            return false;
        }

        double[][] tempGrid = other.getEntries();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (this.entries[i][j] > tempGrid[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isSymmetric() {
        for (int i = 0; i < this.rows; i++) {
            for (int j = i+1; j < this.columns; j++) {
                if (!almostEqual(entries[i][j], entries[j][i])) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isSkewSymmetric() {
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

    private List<Double> getLowerTriangle() {
        if (!this.isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        List<Double> result = new ArrayList<>();
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                if (i >= j) {
                    result.add(entries[i][j]);
                }
            }
        }
        return result;
    }

    private List<Double> getUpperTriangle() {
        if (!this.isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        List<Double> result = new ArrayList<>();
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                if (i <= j) {
                    result.add(entries[i][j]);
                }
            }
        }
        return result;
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
        if (len1 != len2) {
            throw new IllegalArgumentException(String.format("Dot product is undefined for array lengths %d and %d", len1, len2));
        }

        double res = 0;
        for (int i = 0; i < len1; i++) {
            res += (arr1[i] * arr2[i]);
        }
        return res;
    }

    private double[] getColumn(int n) {

        if (!colInRange(n)) {
            throw new IndexOutOfBoundsException(String.format("Column index %d is out of bounds for Matrix of order %s", n, order));
        }
        double[] col = new double[this.rows];
        for (int i = 0; i < this.rows; i++) {
            col[i] = this.getValue(i, n);
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
        if (grid.length == 0) {
            throw new IllegalArgumentException("Matrix grid must be non-empty.");
        }

        int refSize = grid[0].length;
        for (int i = 1; i < grid.length; i++) {
            if (grid[i].length != refSize) {
                return i;
            }
        }
        return -1; //No such row found
    }

    private static int getMismatchedRowIndex(List<List<Double>> grid) {
        if (grid == null) {
            throw new IllegalArgumentException("Matrix grid must be non-null.");
        }

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
        if (other == null || !this.getOrder().equals(other.getOrder())) return false;

        double[][] tempGrid = other.getEntries();
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                if (!almostEqual(this.entries[i][j], tempGrid[i][j])) {
                    return false;
                }
            }
        }
        return true;
    }

    private boolean rowInRange(int rowIndex) {
        return 0 <= rowIndex && rowIndex < this.rows;
    }

    private boolean colInRange(int colIndex) {
        return 0 <= colIndex && colIndex < this.columns;
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

        if (0 > colIdx || colIdx >= grid[0].length) {
            throw new IndexOutOfBoundsException("Column index out of bounds.");
        }

        if (!almostEqual(grid[colIdx][colIdx], 0.0)) {
            return new Pivot(grid[colIdx][colIdx], colIdx, colIdx, false); //Ideal case: pivot is found along diagonal
        }

        for (int rowIdx = colIdx + 1; rowIdx < grid.length; rowIdx++) { //already checked diagonal
            double possibleValue = grid[rowIdx][colIdx];
            if (!almostEqual(possibleValue, 0.0)) {
                return new Pivot(possibleValue, rowIdx, colIdx, true);
                //in this case, row swapping needs to be done until the pivot ends up on the diagonal.
                // track number of row swaps for final det calculation
            }
        }
        return new Pivot(); //valid pivot not found in that column. If a valid pivot is not found, the determiInfinite/undefined valuet is zero
    }

    private void swapGridRow(double[][] grid, int c, int row) {
        double[] temp = grid[c];
        grid[c] = grid[row];
        grid[row] = temp;
    }

    private boolean isDiagonalElement(double value) {
        for (int i = 0; i < this.rows; i++) {
            double diagonalElement = this.entries[i][i];
            if (almostEqual(diagonalElement, value)) {
                return true;
            }
        }
        return false;
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

    private int getArrayHashCode(double[][] grid) {
        int hash = 0;
        for (double[] row : grid) {
            for (double value : row) {
                double normalized = Math.round(value / TOLERANCE) * TOLERANCE;
                hash = (31 * hash) + Double.hashCode(normalized);
            }
        }
        return hash;
    }
    
    private static boolean containsNullRows(double[][] grid) {
        for (double[] row : grid) {
            if (row == null) {
                return true;
            }
        }
        return false;
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
        return 31 * order.hashCode() + getArrayHashCode(this.entries);
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof Matrix)) return false;

        Matrix o = (Matrix)other;
        return this.equalsMatrix(o);
    }
}

