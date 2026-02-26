package main;
import java.util.Arrays;
import java.util.List;

public final class Matrix {

    // ==== INSTANCE VARIABLES AND DECLARATIONS ====

    private static final double TOLERANCE = 1e-6;

    private final int rows;
    private final int columns;
    private final Pair order;
    private final double[][] entries;

    // ==== MATRIX CREATION METHODS ====

    /**
     * Constructs a matrix with the specified number of rows and columns.
     * All entries are initialized to {@code 0.0}.
     *
     * This constructor allocates a new Matrix instance and does not mutate
     * any existing Matrix.
     *
     * @param nrows the number of rows; must be strictly positive.
     * @param ncols the number of columns; must be strictly positive.
     *
     * @throws MatrixException if {@code nrows} <= 0 or {@code ncols} <= 0
     *
     * @return a newly constructed matrix of dimensions {@code nrows} x {@code ncols}
     *
     */
    public Matrix(int nrows, int ncols) {
        if (nrows <= 0 || ncols <= 0) {
            throw MatrixException.illegalDimensions();
        }

        this.rows = nrows;
        this.columns = ncols;
        this.order = new Pair(rows, columns);
        this.entries = new double[rows][columns];
    }


    /**
     * Constructs a matrix from a two-dimensional primitive {@code double} array.
     *
     * The input array is validated to ensure that:
     * 2D array references are not {@code null},
     * no row reference is {@code null},
     * the array contains at least one row,
     * each row contains at least one column,
     * all rows have identical lengths,
     * no element is {@code NaN} or infinite.
     *
     * A deep copy of the provided array is created. Subsequent modification
     * of the input array does not affect this matrix. This constructor does
     * not mutate any external data structure.
     *
     * @param grid a non-null rectangular 2D array of {@code double} values with at least
     *             one row and one column and no {@code NaN} or infinite values
     *
     * @throws IllegalArgumentException if grid is {@code null}, contains a {@code null} row,
     *                                  or contains rows of unequal length
     * @throws MatrixException if {@code grid.length == 0}, any row has length {@code 0},
     *                         or any element is {@code NaN} or infinite
     *
     * @return a newly constructed matrix of dimensions
     *         {@code grid.length} x {@code grid[0].length}
     *
     * Time Complexity: O(m * n), where m = {@code grid.length} and
     * n = {@code grid[0].length}
     */
    public Matrix(double[][] grid) {

        validateGrid(grid);
        this.entries = deepGridCopy(grid);
        this.rows = grid.length;
        this.columns = grid[0].length;
        this.order = new Pair(rows, columns);
    }

    /**
     * Constructs a matrix from a two-dimensional nested {@code List<Double>}.
     *
     * The input structure is validated to ensure:
     * the outer list reference is non-null,
     * no row list reference is null,
     * the outer list contains at least one row,
     * each row contains at least one element,
     * all rows have identical sizes,
     * no element is null,
     * no element is {@code NaN} or infinite.
     *
     * A deep copy of the provided data is created and stored internally.
     * Subsequent modification of the input lists does not affect this matrix.
     * This constructor does not mutate the provided lists.
     *
     * @param grid a non-null rectangular {@code List<List<Double>>}
     *             containing at least one row and one column, with no
     *             null, {@code NaN}, or infinite values
     *
     * @throws IllegalArgumentException if {@code grid} is {@code null},
     *                                  contains a {@code null} row,
     *                                  or contains a {@code null} element
     * @throws MatrixException if {@code grid} is empty,
     *                         any row is empty,
     *                         rows have unequal sizes,
     *                         or any value is {@code NaN} or infinite
     *
     * @return a newly constructed matrix of dimensions
     *         {@code grid.size()} x {@code grid.get(0).size()}
     */
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
        this.columns = grid.get(0).size();
        this.order = new Pair(rows, columns);
        this.entries = new double[rows][columns];

        for (int i = 0; i < rows; i++) {
            List<Double> ithRow = grid.get(i);
            for (int j = 0; j < columns; j++) {
                double value = ithRow.get(j);

                if (!Double.isFinite(value)) {
                    throw new IllegalArgumentException("Infinite/undefined value at coordinates " + new Pair(i,j));
                }
                this.entries[i][j] = value;
            }
        }
    }

    /**
     * Constructs a square matrix with the specified dimension.
     *
     * This constructor initializes a matrix with {@code nrows} rows and
     * {@code nrows} columns. All entries are initialized to {@code 0.0}.
     * Internally delegates to the primary constructor.
     *
     * This constructor creates a new matrix instance and does not mutate
     * any existing matrix.
     *
     * @param nrows the number of rows and columns; must be strictly positive
     *
     * @throws MatrixException if {@code nrows <= 0}
     *
     * @return a newly constructed square matrix of dimensions
     *         nrows × nrows
     */
    public Matrix(int nrows) {
        this(nrows, nrows);
    }

    /**
     * The copy constructor for existing, non-null Matrix instances. It is used for the quick deep copying of existing Matrix
     * instances. Since constructors strictly enforce the creation of valid Matrix instances, we only check for {@code null}
     * references before initializing the copy.
     * To ensure accurate deep copying, we use a deep copy of the existing Matrix instance's {@code entries}
     * (through the {@code toArray()} method) and pass that into the primitive grid constructor.
     *
     * @throws IllegalArgumentException if the existing Matrix is a null reference.
     * @param A the Matrix instance that needs to be deep copied.
     */
    public Matrix(Matrix A) {
        this(getValidMatrix(A));
    }


    /**
     * Creates a matrix from a variable number of primitive {@code double[]} rows.
     *
     * The provided rows are interpreted as the rows of the resulting matrix.
     * The input is validated to ensure that:
     * the outer array reference is non-null,
     * no row reference is {@code null},
     * at least one row is provided,
     * each row has positive length,
     * all rows have identical lengths,
     * no element is {@code NaN} or infinite.
     *
     * A deep copy of the provided rows is created. Subsequent modification
     * of the input arrays does not affect the returned matrix. This method
     * does not mutate the provided arrays.
     *
     * @param rows a non-null sequence of {@code double[]} arrays representing
     *             matrix rows; must form a non-empty rectangular structure
     *             with no {@code NaN} or infinite values.
     *
     * @return a new matrix of dimensions {@code rows.length} × {@code rows[0].length}
     *
     * @throws IllegalArgumentException if {@code rows} is null or
     *                                  contains a null row
     * @throws MatrixException if no rows are provided,
     *                         any row has length 0,
     *                         rows have unequal lengths,
     *                         or any value is {@code NaN} or infinite
     */
    public static Matrix ofRows(double[]... rows) {
        return new Matrix(rows);
    }

    /**
     * Creates a matrix from a variable number of {@code double[]} arrays,
     * interpreted as columns of the resulting matrix.
     *
     * The input must be non-null, contain at least one column, form a
     * rectangular structure with equal non-zero lengths, and contain
     * no {@code NaN} or infinite values.
     *
     * A deep copy of the provided column data is created. The input
     * arrays are not mutated.
     *
     * If {@code k} columns of length {@code m} are provided, the
     * resulting matrix has dimensions {@code m × k}.
     *
     * @param columns non-null column arrays forming a non-empty
     *                rectangular matrix
     *
     * @return a new matrix constructed from the given columns
     *
     * @throws IllegalArgumentException if {@code columns} is null
     *                                  or contains a null column reference
     * @throws MatrixException if no columns are provided,
     *                         column lengths differ or are zero,
     *                         or any value is {@code NaN} or infinite
     */
    public static Matrix ofColumns(double[]... columns) {
        return Matrix.ofRows(columns).transpose();
    }


    /**
     * Creates a matrix of specified dimensions in which every entry is
     * initialized to the scalar {@code k}.
     *
     * The dimensions must be strictly positive. The scalar must be finite.
     * If {@code |k| < TOLERANCE} (in this case, {@code 1e-6}), the resulting matrix is
     * treated as a zero matrix.
     *
     * This method returns a new matrix instance and does not mutate
     * any existing matrix.
     *
     * @param nrows the number of rows; must be strictly positive
     * @param ncols the number of columns; must be strictly positive
     * @param k the finite scalar value assigned to each entry
     *
     * @return a new {@code nrows × ncols} matrix with all entries equal to {@code k},
     *         or the zero matrix if {@code |k| < TOLERANCE}
     *
     * @throws MatrixException if {@code nrows <= 0}, {@code ncols <= 0},
     *                         or {@code k} is {@code NaN} or infinite
     */
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


    /**
     * Mutates this matrix by setting all entries to {@code 0.0}.
     *
     * The dimensions of the matrix remain unchanged.
     * This operation does not allocate a new matrix.
     *
     * Time Complexity: O(m * n)
     */
    public void zeroMatrixInPlace() {
        fillInPlace(0.0);
    }


    /**
     * Creates a square scalar matrix of dimension {@code nrows}.
     *
     * The resulting matrix has {@code k} on its main diagonal and
     * {@code 0.0} on all off-diagonal entries. The dimension must be
     * strictly positive and {@code k} must be finite.
     *
     * If {@code |k| < TOLERANCE} (in this case, {@code 1e-6}), the resulting matrix
     * is treated as the zero matrix.
     *
     * This method returns a new matrix instance and does not mutate
     * any existing matrix.
     *
     * @param nrows the number of rows and columns; must be strictly positive
     * @param k the finite scalar value assigned to each diagonal entry
     *
     * @return a new {@code nrows × nrows} scalar matrix, or the zero matrix
     *         if {@code |k| < TOLERANCE}
     *
     * @throws MatrixException if {@code nrows <= 0} or {@code k} is {@code NaN} or infinite
     */
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


    /**
     * Creates an {@code nrows × nrows} identity matrix.
     *
     * The resulting matrix has {@code 1.0} on its main diagonal and
     * {@code 0.0} on all off-diagonal entries. The dimension must be
     * strictly positive.
     *
     * This method returns a new matrix instance and does not mutate
     * any existing matrix.
     *
     * @param nrows the number of rows and columns; must be strictly positive
     *
     * @return a new {@code nrows × nrows} identity matrix
     *
     * @throws MatrixException if {@code nrows <= 0}
     */
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


    /**
     * Returns the entry at the specified row and column.
     *
     * Indices follow 0-based indexing.
     * This method does not mutate the matrix.
     *
     * @param rowIndex the row index; must satisfy
     *                 {@code 0 <= rowIndex < rows}
     * @param colIndex the column index; must satisfy
     *                 {@code 0 <= colIndex < columns}
     *
     * @return the value at position {@code (rowIndex, colIndex)}
     *
     * @throws IndexOutOfBoundsException if {@code rowIndex} or
     *                                   {@code colIndex} is outside
     *                                   the valid index range
     */
    public double getEntry(int rowIndex, int colIndex) {
        return getValue(rowIndex, colIndex);
    }


    /**
     * Sets the entry at position {@code (r, c)} to the specified value.
     *
     * The value must be finite. Indices follow 0-based indexing and must
     * refer to a valid position in the matrix.
     *
     * This operation mutates the matrix in place and modifies only the
     * specified entry.
     *
     * @param value the finite value to assign to position {@code (r, c)}
     * @param r the row index; must satisfy {@code 0 <= r < rows}
     * @param c the column index; must satisfy {@code 0 <= c < columns}
     *
     * @throws IllegalArgumentException if {@code value} is {@code NaN}
     *                                  or infinite
     * @throws IndexOutOfBoundsException if {@code r} or {@code c}
     *                                   is outside the valid index range
     */
    public void setEntry(double value, int r, int c) {
        if (!isInBounds(r,c)) {
            throw new IndexOutOfBoundsException(String.format("Cell coordinates (%d, %d) out of bounds for Matrix of order %s", r, c, order));
        }

        if (!Double.isFinite(value)) {
            throw new IllegalArgumentException("Given value cannot be used in matrix since it is Not a Number (Infinite/undefined value)");
        }

        this.entries[r][c] = value;
    }


    /**
     * Replaces the row at {@code rowIndex} with the specified values.
     *
     * The provided array must be non-null, have length equal to the
     * number of columns in the matrix, and contain only finite values.
     *
     * This operation mutates the matrix in place and does not allocate
     * additional matrix storage.
     *
     * @param row the non-null array of finite values whose length must
     *            equal {@code columns}
     * @param rowIndex the 0-based index of the row to replace; must satisfy
     *                 {@code 0 <= rowIndex < rows}
     *
     * @throws IllegalArgumentException if {@code row} is null or contains
     *                                  {@code NaN} or infinite values
     * @throws MatrixException if {@code row.length != getColumnCount()}
     * @throws IndexOutOfBoundsException if {@code rowIndex} is outside
     *                                   the valid index range
     */
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

    /**
     * Replaces the column at {@code colIndex} with the specified values.
     *
     * The provided array must be non-null, have length equal to the
     * number of rows in the matrix, and contain only finite values.
     *
     * This operation mutates the matrix in place and does not allocate
     * additional matrix storage.
     *
     * @param col the non-null array of finite values whose length must
     *            equal {@code rows}
     * @param colIndex the 0-based index of the column to replace; must satisfy
     *                 {@code 0 <= colIndex < rows}
     *
     * @throws IllegalArgumentException if {@code col} is null or contains
     *                                  {@code NaN} or infinite values
     * @throws MatrixException if {@code col.length != rows}
     * @throws IndexOutOfBoundsException if {@code colIndex} is outside
     *                                   the valid index range
     */
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


    /**
     * Scales each entry of this matrix by the finite real scalar {@code k}.
     *
     * This operation mutates {@code this} in-place and does not allocate
     * a new matrix. The dimensions remain unchanged.
     *
     * @param k the scalar factor applied to each entry
     * @throws MatrixException if {@code k} is infinite or {@code NaN}
     */
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

    /**
     * Adds the entries of {@code other} to this matrix element-wise.
     *
     * Both matrices must be non-null and have identical dimensions.
     *
     * This operation mutates this matrix in place. The {@code other}
     * matrix remains unchanged.
     *
     * @param other the matrix to be added; must be non-null and have
     *              the same dimensions as this matrix
     *
     * @throws IllegalArgumentException if {@code other} is null
     * @throws MatrixException if the dimensions of the matrices differ
     *
     * Time Complexity: O(m × n), where m = number of rows and
     * n = number of columns
     */
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

    /**
     * Subtracts the entries of {@code other} from this matrix element-wise.
     *
     * Both matrices must be non-null and have identical dimensions.
     *
     * This operation mutates this matrix in place. The {@code other}
     * matrix remains unchanged.
     *
     * @param other the matrix to be subtracted; must be non-null and have
     *              the same dimensions as this matrix
     *
     * @throws IllegalArgumentException if {@code other} is null
     * @throws MatrixException if the dimensions of the matrices differ
     *
     * Time Complexity: O(m × n), where m = number of rows and
     * n = number of columns
     */
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

    /**
     * Transposes this matrix in place by reflecting it across its main diagonal.
     *
     * The matrix must be square. Each pair of symmetric off-diagonal
     * elements {@code (i, j)} and {@code (j, i)} is swapped exactly once.
     *
     * This operation mutates the matrix in place and does not allocate
     * additional matrix storage.
     *
     * @throws MatrixException if the matrix is not square
     *
     * Time Complexity: O(n^2), where n = number of rows (and columns)
     */
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

    public Matrix symmetricPart() {
        if (!isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        Matrix result = this.add(this.transpose());
        result.multiplyByScalarInPlace(0.5);
        return result;
    }

    public Matrix skewSymmetricPart() {
        if (!isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        Matrix result = this.subtract(this.transpose());
        result.multiplyByScalarInPlace(0.5);
        return result;
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

    private static boolean containsNullRows(List<List<Double>> grid) {
        for (List<Double> row : grid) {
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

