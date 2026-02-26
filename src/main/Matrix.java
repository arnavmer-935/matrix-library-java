package main;
import java.util.Arrays;
import java.util.List;

/**
 * A dense, row-major implementation of a real-valued matrix
 * supporting core linear algebra operations and structural queries.
 *
 * <p>This class provides:
 * <ul>
 *   <li>Arithmetic operations (addition, subtraction, scalar multiplication, multiplication)</li>
 *   <li>Linear algebra operations (determinant, inverse, transpose)</li>
 *   <li>Matrix decompositions (symmetric and skew-symmetric parts)</li>
 *   <li>Structural and property queries (e.g., square, triangular, symmetric)</li>
 * </ul>
 *
 * <p><strong>Storage Model:</strong><br>
 * Matrices are stored in row-major order using a {@code double[][]}
 * backing array.
 *
 * <p><strong>Mutability Design:</strong><br>
 * Operations are provided in both in-place and out-of-place variants.
 * In-place methods mutate the current instance, while non in-place
 * methods return new {@code Matrix} instances and leave the original
 * unchanged.
 *
 * <p><strong>Numerical Policy:</strong><br>
 * Floating-point comparisons in structural queries and pivot
 * detection use a configured numerical tolerance to account for
 * rounding behavior. However, {@link #equals(Object)} performs
 * strict element-wise comparison using exact {@code double}
 * semantics to preserve the {@code equals}/{@code hashCode}
 * contract.
 *
 * <p><strong>Validation Guarantees:</strong><br>
 * All dimension-sensitive operations validate input sizes.
 * Linear algebra operations that require specific matrix configurations
 * will throw a {@code MatrixException} if the precondition is not met.
 * Any invalid arguments or user inputs (such as uneven row lengths,
 * negative dimensions, etc.) will cause the program to throw an
 * {@code IllegalArgumentException}.
 * For indexing inaccuracies, the program throws a standard
 * {@code IndexOutOfBoundsException}
 *
 * <p>This class is designed for correctness, clarity, and explicit
 * numerical behavior rather than high-performance optimized
 * computation.
 */
public final class Matrix {

    /**
     *
     * All Matrix instances guarantee the following information about its attributes:
     * - rows > 0
     * - columns > 0
     * - entries != null
     * - entries.length == rows
     * - For every r in [0, rows):
     *       entries[r] != null
     *       entries[r].length == columns
     *
     * - entries is deep-owned by this instance.
     *   No external references to internal row arrays exist.
     *
     * - Equality (equals/hashCode) uses strict element-wise
     *   comparison and does NOT use tolerance.
     *
     * - Numerical queries (e.g., structural checks) use the
     *   configured equality tolerance.
     *
     * - Pivot detection in row-reduction and inverse uses a
     *   stricter pivot tolerance to handle near-singular matrices.
     */
    private static final double TOLERANCE = 1e-6;

    private final int rows;
    private final int columns;
    private final Pair order;
    private final double[][] entries;

    // ==== MATRIX CREATION METHODS ====

    /**
     * Constructs a {@code Matrix} with the specified number of rows and columns.
     * All entries are initialized to {@code 0.0}.
     *
     * <p>Both {@code nrows} and {@code ncols} must be strictly positive.
     * The resulting matrix has dimension {@code nrows × ncols}.
     *
     * <p>This constructor allocates a new {@code Matrix} instance.
     * No existing matrix is mutated.
     *
     * <p><strong>Time Complexity:</strong> O(nrows · ncols)<br>
     * <strong>Space Complexity:</strong> O(nrows · ncols)
     *
     * @param nrows the number of rows; must be greater than {@code 0}
     * @param ncols the number of columns; must be greater than {@code 0}
     *
     * @throws MatrixException if {@code nrows <= 0} or {@code ncols <= 0}
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
     * Constructs a {@code Matrix} from a two-dimensional primitive {@code double} array.
     *
     * <p>The input array must be non-{@code null}, contain at least one row,
     * and each row must be non-{@code null} with at least one column.
     * All rows must have identical length, forming a rectangular structure.
     * Every element must be finite (no {@code NaN} or infinite values).
     *
     * <p>A deep copy of {@code grid} is created. Subsequent modifications to
     * the input array do not affect this matrix. No external data structure
     * is mutated.
     *
     * <p>The resulting matrix has dimension
     * {@code grid.length × grid[0].length}.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(m · n)<br>
     * where {@code m = grid.length} and
     * {@code n = grid[0].length}.
     *
     * @param grid a rectangular 2D array of finite {@code double} values
     *
     * @throws IllegalArgumentException if {@code grid} is {@code null},
     *         contains a {@code null} row, or contains rows of unequal length
     * @throws MatrixException if {@code grid.length == 0},
     *         any row has length {@code 0}, or any element is
     *         {@code NaN} or infinite
     */
    public Matrix(double[][] grid) {

        validateGrid(grid);
        this.entries = deepGridCopy(grid);
        this.rows = grid.length;
        this.columns = grid[0].length;
        this.order = new Pair(rows, columns);
    }

    /**
     * Constructs a {@code Matrix} from a two-dimensional nested {@code List<Double>}.
     *
     * <p>The outer list must be non-{@code null} and contain at least one row.
     * Each row must be non-{@code null}, contain at least one element,
     * and all rows must have identical sizes, forming a rectangular structure.
     * Elements must be non-{@code null} and finite (no {@code NaN} or infinite values).
     *
     * <p>A deep copy of the provided data is created and stored internally.
     * Subsequent modification of the input lists does not affect this matrix.
     * No provided list is mutated.
     *
     * <p>The resulting matrix has dimension
     * {@code grid.size() × grid.get(0).size()}.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(m · n)<br>
     * where {@code m = grid.size()} and
     * {@code n = grid.get(0).size()}.
     *
     * @param grid a rectangular {@code List<List<Double>>} containing
     *             finite, non-{@code null} values
     *
     * @throws IllegalArgumentException if {@code grid} is {@code null},
     *         contains a {@code null} row, or contains a {@code null} element
     * @throws MatrixException if {@code grid} is empty,
     *         any row is empty, rows have unequal sizes,
     *         or any value is {@code NaN} or infinite
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
     * Constructs a square {@code Matrix} of dimension {@code nrows × nrows}.
     * All entries are initialized to {@code 0.0}.
     *
     * <p>The dimension {@code nrows} must be strictly positive.
     * This constructor delegates to the primary rectangular constructor.
     *
     * <p>A new {@code Matrix} instance is allocated.
     * No existing matrix is mutated.
     *
     * <p><strong>Time Complexity:</strong> O(nrows²)<br>
     * <strong>Space Complexity:</strong> O(nrows²)
     *
     * @param nrows the number of rows and columns; must be greater than {@code 0}
     *
     * @throws MatrixException if {@code nrows <= 0}
     */
    public Matrix(int nrows) {
        this(nrows, nrows);
    }

    /**
     * Constructs a new {@code Matrix} as a deep copy of an existing instance.
     *
     * <p>The provided matrix {@code A} must be non-{@code null}.
     * Since all {@code Matrix} constructors enforce structural and numerical
     * validity, no additional validation beyond a {@code null} check is required.
     *
     * <p>A deep copy of {@code A}'s internal entries is created via
     * {@code A.toArray()} and delegated to the primitive-array constructor.
     * The resulting matrix is structurally independent of {@code A}.
     *
     * <p>This constructor allocates a new {@code Matrix} instance.
     * The provided matrix is not mutated.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(m · n)<br>
     * where {@code m} and {@code n} are the dimensions of {@code A}.
     *
     * @param A the matrix to be deep copied
     *
     * @throws IllegalArgumentException if {@code A} is {@code null}
     */
    public Matrix(Matrix A) {
        this(getValidMatrix(A));
    }


    /**
     * Creates a {@code Matrix} from a variable number of primitive
     * {@code double[]} row arrays.
     *
     * <p>The provided arrays are interpreted as the rows of the resulting matrix.
     * The outer array must be non-{@code null} and contain at least one row.
     * Each row must be non-{@code null}, have positive length, and all rows
     * must have identical length, forming a rectangular structure.
     * Every element must be finite (no {@code NaN} or infinite values).
     *
     * <p>A deep copy of the provided row arrays is created.
     * Subsequent modification of the input arrays does not affect the
     * constructed matrix. No provided array is mutated.
     *
     * <p>The resulting matrix has dimension
     * {@code rows.length × rows[0].length}.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(m · n)<br>
     * where {@code m = rows.length} and
     * {@code n = rows[0].length}.
     *
     * @param rows a sequence of {@code double[]} arrays representing
     *             matrix rows
     *
     * @throws IllegalArgumentException if {@code rows} is {@code null}
     *         or contains a {@code null} row
     * @throws MatrixException if no rows are provided,
     *         any row has length {@code 0},
     *         rows have unequal lengths,
     *         or any value is {@code NaN} or infinite
     */
    public static Matrix ofRows(double[]... rows) {
        return new Matrix(rows);
    }

    /**
     * Creates a {@code Matrix} from a variable number of {@code double[]}
     * arrays interpreted as columns.
     *
     * <p>The outer array must be non-{@code null} and contain at least one column.
     * Each column must be non-{@code null}, have positive length, and all columns
     * must have identical length, forming a rectangular structure.
     * Every element must be finite (no {@code NaN} or infinite values).
     *
     * <p>If {@code k} columns of length {@code m} are provided,
     * the resulting matrix has dimension {@code m × k}.
     *
     * <p>A deep copy of the provided column arrays is created.
     * The input arrays are not mutated.
     *
     * <p><strong>Time Complexity:</strong> O(m · k)<br>
     * <strong>Space Complexity:</strong> O(m · k)<br>
     * where {@code m} is the column length and {@code k} is the number of columns.
     *
     * @param columns column arrays forming a rectangular matrix
     *
     * @throws IllegalArgumentException if {@code columns} is {@code null}
     *         or contains a {@code null} column reference
     * @throws MatrixException if no columns are provided,
     *         column lengths are zero or unequal,
     *         or any value is {@code NaN} or infinite
     */
    public static Matrix ofColumns(double[]... columns) {
        return Matrix.ofRows(columns).transpose();
    }


    /**
     * Creates a {@code Matrix} of dimension {@code nrows × ncols}
     * in which every entry is initialized to the scalar {@code k}.
     *
     * <p>Both {@code nrows} and {@code ncols} must be strictly positive.
     * The scalar {@code k} must be finite.
     *
     * <p>If {@code |k| <= TOLERANCE} (e.g., {@code 1e-6}),
     * the resulting matrix is treated as a zero matrix.
     *
     * <p>A new {@code Matrix} instance is allocated.
     * No existing matrix is mutated.
     *
     * <p><strong>Time Complexity:</strong> O(nrows · ncols)<br>
     * <strong>Space Complexity:</strong> O(nrows · ncols)
     *
     * @param nrows the number of rows; must be greater than {@code 0}
     * @param ncols the number of columns; must be greater than {@code 0}
     * @param k the finite scalar assigned to each entry
     *
     * @return a new {@code nrows × ncols} matrix with all entries equal to {@code k},
     *         or the zero matrix if {@code |k| <= TOLERANCE}
     *
     * @throws MatrixException if {@code nrows <= 0},
     *         {@code ncols <= 0}, or {@code k} is {@code NaN} or infinite
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
     * Sets every entry of {@code this} matrix to {@code 0.0}.
     *
     * <p>This operation <strong>mutates</strong> {@code this}.
     * The matrix dimensions remain unchanged.
     * No new {@code Matrix} instance is allocated.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(1)<br>
     * where {@code m} and {@code n} are the dimensions of {@code this}.
     */
    public void zeroMatrixInPlace() {
        fillInPlace(0.0);
    }

    /**
     * Creates a square scalar {@code Matrix} of dimension {@code nrows × nrows}.
     *
     * <p>The resulting matrix has {@code k} on its main diagonal and
     * {@code 0.0} on all off-diagonal entries.
     * The dimension {@code nrows} must be strictly positive,
     * and {@code k} must be finite.
     *
     * <p>If {@code |k| <= TOLERANCE} (e.g., {@code 1e-6}),
     * the resulting matrix is treated as the zero matrix.
     *
     * <p>This method allocates and returns a new {@code Matrix} instance.
     * No existing matrix is mutated.
     *
     * <p><strong>Time Complexity:</strong> O(nrows²)<br>
     * <strong>Space Complexity:</strong> O(nrows²)
     *
     * @param nrows the number of rows and columns; must be greater than {@code 0}
     * @param k the finite scalar assigned to each diagonal entry
     *
     * @return a new {@code nrows × nrows} scalar matrix,
     *         or the zero matrix if {@code |k| <= TOLERANCE}
     *
     * @throws MatrixException if {@code nrows <= 0}
     *         or {@code k} is {@code NaN} or infinite
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
     * Creates an identity {@code Matrix} of dimension {@code nrows × nrows}.
     *
     * <p>The resulting matrix has {@code 1.0} on its main diagonal
     * and {@code 0.0} on all off-diagonal entries.
     * The dimension {@code nrows} must be strictly positive.
     *
     * <p>This method allocates and returns a new {@code Matrix} instance.
     * No existing matrix is mutated.
     *
     * <p><strong>Time Complexity:</strong> O(nrows²)<br>
     * <strong>Space Complexity:</strong> O(nrows²)
     *
     * @param nrows the number of rows and columns; must be greater than {@code 0}
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
     * Returns the entry at position {@code (rowIndex, colIndex)}.
     *
     * <p>Indices use 0-based indexing.
     * Valid indices satisfy
     * {@code 0 <= rowIndex < rows} and
     * {@code 0 <= colIndex < columns}.
     *
     * <p>This method does <strong>not</strong> mutate {@code this}.
     * @param rowIndex the row index
     * @param colIndex the column index
     *
     * @return the value at the specified position
     *
     * @throws IndexOutOfBoundsException if {@code rowIndex}
     *         or {@code colIndex} is outside the valid index range
     */
    public double getEntry(int rowIndex, int colIndex) {
        return getValue(rowIndex, colIndex);
    }

    /**
     * Sets the entry at position {@code (r, c)} to the specified value.
     *
     * <p>Indices use 0-based indexing and must satisfy
     * {@code 0 <= r < rows} and {@code 0 <= c < columns}.
     * The provided {@code value} must be finite.
     *
     * <p>This operation <strong>mutates</strong> {@code this}
     * and modifies only the specified entry.
     *
     * @param value the finite value to assign
     * @param r the row index
     * @param c the column index
     *
     * @throws IllegalArgumentException if {@code value} is {@code NaN}
     *         or infinite
     * @throws IndexOutOfBoundsException if {@code r} or {@code c}
     *         is outside the valid index range
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
     * <p>The provided {@code row} must be non-{@code null}, have length equal to
     * the number of columns in {@code this}, and contain only finite values.
     * The index {@code rowIndex} must satisfy
     * {@code 0 <= rowIndex < rows}.
     *
     * <p>This operation <strong>mutates</strong> {@code this} in place.
     * No additional matrix storage is allocated.
     *
     * @param row the replacement row values
     * @param rowIndex the 0-based index of the row to replace
     *
     * @throws IllegalArgumentException if {@code row} is {@code null}
     *         or contains {@code NaN} or infinite values
     * @throws MatrixException if {@code row.length != getColumnCount()}
     * @throws IndexOutOfBoundsException if {@code rowIndex}
     *         is outside the valid index range
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
     * <p>The provided {@code col} must be non-{@code null}, have length equal to
     * the number of rows in {@code this}, and contain only finite values.
     * The index {@code colIndex} must satisfy
     * {@code 0 <= colIndex < columns}.
     *
     * <p>This operation <strong>mutates</strong> {@code this} in place.
     * No additional matrix storage is allocated.
     *
     * @param col the replacement column values
     * @param colIndex the 0-based index of the column to replace
     *
     * @throws IllegalArgumentException if {@code col} is {@code null}
     *         or contains {@code NaN} or infinite values
     * @throws MatrixException if {@code col.length != getRowCount()}
     * @throws IndexOutOfBoundsException if {@code colIndex}
     *         is outside the valid index range
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
     * Scales every entry of {@code this} matrix by the finite scalar {@code k}.
     *
     * <p>The scalar {@code k} must be finite.
     *
     * <p>This operation <strong>mutates</strong> {@code this} in place.
     * The matrix dimensions remain unchanged.
     * No new {@code Matrix} instance is allocated.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(1)<br>
     * where {@code m} and {@code n} are the dimensions of {@code this}.
     *
     * @param k the scalar factor applied to each entry
     *
     * @throws MatrixException if {@code k} is {@code NaN} or infinite
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
     * Adds the entries of {@code other} to {@code this} matrix element-wise.
     *
     * <p>The provided matrix {@code other} must be non-{@code null}
     * and have identical dimensions to {@code this}.
     *
     * <p>This operation <strong>mutates</strong> {@code this} in place.
     * The {@code other} matrix remains unchanged.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(1)<br>
     * where {@code m} and {@code n} are the dimensions of {@code this}.
     *
     * @param other the matrix to be added
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}
     * @throws MatrixException if the dimensions of {@code this}
     *         and {@code other} differ
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
     * Subtracts the entries of {@code other} from {@code this} matrix element-wise.
     *
     * <p>The provided matrix {@code other} must be non-{@code null}
     * and have identical dimensions to {@code this}.
     *
     * <p>This operation <strong>mutates</strong> {@code this} in place.
     * The {@code other} matrix remains unchanged.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(1)<br>
     * where {@code m} and {@code n} are the dimensions of {@code this}.
     *
     * @param other the matrix to be subtracted
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}
     * @throws MatrixException if the dimensions of {@code this}
     *         and {@code other} differ
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
     * Transposes {@code this} matrix in place by reflecting it across
     * its main diagonal.
     *
     * <p>The matrix must be square.
     * Each pair of symmetric off-diagonal entries
     * {@code (i, j)} and {@code (j, i)} is swapped exactly once.
     *
     * <p>This operation <strong>mutates</strong> {@code this}.
     * No additional matrix storage is allocated.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(1)<br>
     * where {@code n} is the dimension of the matrix.
     *
     * @throws MatrixException if {@code this} is not square
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
    /**
     * Scales each entry of this matrix by the finite real scalar {@code k},
     * and returns a new Matrix instance consisting of those scaled entries.
     *
     * This operation does not mutate {@code this} in-place.
     * It allocates a new matrix. The dimensions remain unchanged.
     *
     * @param k the scalar factor applied to each entry
     * @throws MatrixException if {@code k} is infinite or {@code NaN}
     *
     * @return a new Matrix containing entries scaled by {@code k}.
     * Returns a deep copy of the original Matrix if {@code k} is {@code 1.0}.
     * Returns a new instance of a zero Matrix if {@code |k| <= TOLERANCE}.
     */
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

    /**
     * Returns a new {@code Matrix} containing the element-wise sum of
     * {@code this} and {@code other}.
     *
     * <p>The provided matrix {@code other} must be non-{@code null}
     * and have identical dimensions to {@code this}.
     *
     * <p>This operation does <strong>not</strong> mutate {@code this}.
     * The {@code other} matrix also remains unchanged.
     * A new {@code Matrix} instance is allocated to store the result.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(m · n)<br>
     * where {@code m} and {@code n} are the dimensions of {@code this}.
     *
     * @param other the matrix to add
     *
     * @return a new {@code Matrix} containing the element-wise sums
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}
     * @throws MatrixException if the dimensions of {@code this}
     *         and {@code other} differ
     */
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


    /**
     * Returns a new {@code Matrix} containing the element-wise difference
     * of {@code this} and {@code other}.
     *
     * <p>The provided matrix {@code other} must be non-{@code null}
     * and have identical dimensions to {@code this}.
     *
     * <p>This operation does <strong>not</strong> mutate {@code this}.
     * The {@code other} matrix also remains unchanged.
     * A new {@code Matrix} instance is allocated to store the result.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(m · n)<br>
     * where {@code m} and {@code n} are the dimensions of {@code this}.
     *
     * @param other the matrix to subtract
     *
     * @return a new {@code Matrix} containing the element-wise differences
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}
     * @throws MatrixException if the dimensions of {@code this}
     *         and {@code other} differ
     */
    public Matrix subtract(Matrix other) {
        if (other == null) {
            throw new IllegalArgumentException("Matrix operand must be non-null.");
        }
        return add(other.multiplyByScalar(-1));
    }

    /**
     * Returns a new {@code Matrix} representing the transpose of {@code this}.
     *
     * <p>If {@code this} has dimension {@code m × n}, the resulting matrix
     * has dimension {@code n × m}. Each entry satisfies
     * {@code result[j][i] = this[i][j]}.
     * The matrix need not be square.
     *
     * <p>This operation does <strong>not</strong> mutate {@code this}.
     * A new {@code Matrix} instance is allocated to store the transposed entries.
     *
     * <p><strong>Time Complexity:</strong> O(m · n)<br>
     * <strong>Space Complexity:</strong> O(m · n)<br>
     * where {@code m} and {@code n} are the dimensions of {@code this}.
     *
     * @return a new {@code Matrix} representing the transpose of {@code this}
     */
    public Matrix transpose() {
        double[][] transposedGrid = new double[this.columns][this.rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                transposedGrid[j][i] = this.entries[i][j];
            }
        }
        return new Matrix(transposedGrid);
    }


    /**
     * Computes the matrix product of {@code this} and {@code other}
     * and returns a new {@code Matrix} containing the result.
     *
     * <p>The product is defined only if the number of columns in
     * {@code this} equals the number of rows in {@code other}.
     *
     * <p>Multiplication is performed using the standard dot product:
     * each entry {@code (i, j)} in the result is computed as the
     * dot product of row {@code i} of {@code this} and column {@code j}
     * of {@code other}.
     *
     * <p>This operation does <strong>not</strong> mutate {@code this}.
     * A new {@code Matrix} instance is allocated to store the result.
     *
     * <p>If {@code this} has dimension {@code m × n} and
     * {@code other} has dimension {@code n × p},
     * the resulting matrix has dimension {@code m × p}.
     *
     * <p><strong>Time Complexity:</strong> O(m · n · p)<br>
     * <strong>Space Complexity:</strong> O(m · p)
     *
     * <p>If both matrices are {@code n × n}:
     * <br><strong>Time Complexity:</strong> O(n³)
     * <br><strong>Space Complexity:</strong> O(n²)
     *
     * @param other the matrix to multiply with {@code this}
     *
     * @return a new {@code Matrix} representing the product
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}
     * @throws MatrixException if the number of columns in {@code this}
     *         is not equal to the number of rows in {@code other}
     */
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

    /**
     * Returns a deep copy of the internal 2D array representation
     * of this {@code Matrix}.
     *
     * <p>The returned array preserves the row-major ordering of entries
     * but is fully independent of this instance. Modifying the returned
     * array will <strong>not</strong> affect this {@code Matrix}.
     *
     * <p>This method performs a deep copy of all entries to preserve
     * encapsulation and prevent external mutation of internal state.
     *
     * <p><strong>Time Complexity:</strong> O(mn)<br>
     * <strong>Space Complexity:</strong> O(mn)
     *
     * @return a new {@code double[][]} containing the same values as
     *         this {@code Matrix}
     */
    public double[][] toArray() {
        return deepGridCopy(this.entries);
    }

    /**
     * Computes and returns the symmetric part of this {@code Matrix}.
     *
     * <p>The symmetric part is defined as:
     * <pre>
     *     (A + Aᵀ) / 2
     * </pre>
     * where {@code Aᵀ} denotes the transpose of this matrix.
     *
     * <p>This operation does <strong>not</strong> mutate {@code this}.
     * A new {@code Matrix} instance is allocated to store the result.
     *
     * <p>This method is only defined for square matrices.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(n²)
     *
     * @return a new {@code Matrix} representing the symmetric part
     *         of {@code this}
     *
     * @throws MatrixException if this matrix is not square
     */
    public Matrix symmetricPart() {
        if (!isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        Matrix result = this.add(this.transpose());
        result.multiplyByScalarInPlace(0.5);
        return result;
    }

    /**
     * Computes and returns the skew-symmetric part of this {@code Matrix}.
     *
     * <p>The skew-symmetric part is defined as:
     * <pre>
     *     (A - Aᵀ) / 2
     * </pre>
     * where {@code Aᵀ} denotes the transpose of this matrix.
     *
     * <p>This operation does <strong>not</strong> mutate {@code this}.
     * A new {@code Matrix} instance is allocated to store the result.
     *
     * <p>This method is only defined for square matrices.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(n²)
     *
     * @return a new {@code Matrix} representing the skew-symmetric part
     *         of {@code this}
     *
     * @throws MatrixException if this matrix is not square
     */
    public Matrix skewSymmetricPart() {
        if (!isSquareMatrix()) {
            throw MatrixException.requireSquareMatrix();
        }

        Matrix result = this.subtract(this.transpose());
        result.multiplyByScalarInPlace(0.5);
        return result;
    }

    // ==== NUMERICAL METHODS ====

    /**
     * Computes and returns the inverse of this {@code Matrix}
     * using Gauss–Jordan elimination with partial pivoting.
     *
     * <p>The algorithm augments this matrix with the identity matrix
     * and performs row operations until the left side reduces to the
     * identity matrix. The transformed right side then becomes the inverse.
     *
     * <p>Pivot selection uses a numerical tolerance of <strong>{@code 1e-12}</strong> to determine whether
     * a pivot element is sufficiently non-zero. If no valid pivot can be
     * found within this tolerance, the matrix is treated as singular.
     *
     * <p>It should be noted that Pivot Tolerance is different from the actual tolerance
     * of {@code 1e-6} to ensure that leading entries less than {@code 1e-6} can be detected during row reduction.</p>
     *
     * <p>This method does <strong>not</strong> mutate {@code this}.
     * A new {@code Matrix} instance is allocated to store the result.
     *
     * <p>This operation is only defined for square, non-singular matrices.
     *
     * <p><strong>Time Complexity:</strong> O(n³)<br>
     * <strong>Space Complexity:</strong> O(n²)
     *
     * @return a new {@code Matrix} representing the inverse of {@code this}
     *
     * @throws MatrixException if this matrix is not square
     * @throws MatrixException if this matrix is singular (including
     *         cases where pivots fall below the numerical tolerance)
     */
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

    /**
     * Computes and returns the determinant of this {@code Matrix}.
     *
     * <p>This method is only defined for square matrices.
     *
     * <p>For 1 × 1 and 2 × 2 matrices, the determinant is computed
     * directly using closed-form expressions.
     *
     * <p>For larger matrices, the determinant is computed using
     * Gaussian elimination with partial pivoting. The matrix is
     * transformed into an upper triangular form, and the determinant
     * is obtained as the product of the diagonal entries, adjusted
     * for row swaps.
     *
     * <p>Pivot selection uses a numerical tolerance to determine
     * whether a pivot element is sufficiently non-zero. If no valid
     * pivot can be found within this tolerance, the determinant is
     * treated as zero.
     *
     * <p>This operation does <strong>not</strong> mutate {@code this}.
     *
     * <p><strong>Time Complexity:</strong> O(n³)<br>
     * <strong>Space Complexity:</strong> O(n²)
     *
     * @return the determinant of this {@code Matrix}
     *
     * @throws MatrixException if this matrix is not square
     */
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
    /**
     * Returns {@code true} if this {@code Matrix} is square.
     *
     * <p>A matrix is square if its number of rows is equal to
     * its number of columns.
     *
     * <p><strong>Time Complexity:</strong> O(1)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if {@code this.rows == this.columns};
     *         {@code false} otherwise
     */
    public boolean isSquareMatrix() {
        return this.rows == this.columns;
    }

    /**
     * Returns {@code true} if this {@code Matrix} represents
     * the identity matrix.
     *
     * <p>A matrix is considered an identity matrix if:
     * <ul>
     *   <li>It is square.</li>
     *   <li>All diagonal entries are approximately equal to 1.0.</li>
     *   <li>All off-diagonal entries are approximately equal to 0.0.</li>
     * </ul>
     *
     * <p>Comparisons are performed using the class's numerical tolerance of
     * {@code 1e-6} to account for floating-point rounding behavior.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if this matrix is approximately an
     *         identity matrix within the configured tolerance;
     *         {@code false} otherwise
     */
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

    /**
     * Returns {@code true} if this {@code Matrix} is upper triangular.
     *
     * <p>A matrix is considered upper triangular if:
     * <ul>
     *   <li>It is square.</li>
     *   <li>All entries strictly below the main diagonal are
     *       approximately equal to 0.0.</li>
     * </ul>
     *
     * <p>Comparisons are performed using the class's configured
     * numerical tolerance to account for floating-point rounding behavior.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if this matrix is approximately upper
     *         triangular; {@code false} otherwise
     */
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

    /**
     * Returns {@code true} if this {@code Matrix} is lower triangular.
     *
     * <p>A matrix is considered lower triangular if:
     * <ul>
     *   <li>It is square.</li>
     *   <li>All entries strictly above the main diagonal are
     *       approximately equal to 0.0.</li>
     * </ul>
     *
     * <p>Comparisons are performed using the class's configured
     * numerical tolerance to account for floating-point rounding behavior.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if this matrix is approximately lower
     *         triangular; {@code false} otherwise
     */
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

    /**
     * Returns {@code true} if all entries of this {@code Matrix}
     * are approximately equal to 0.0.
     *
     * <p>Comparisons are performed using the class's configured
     * numerical tolerance to account for floating-point rounding behavior.
     *
     * <p><strong>Time Complexity:</strong> O(mn)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if all entries are approximately zero;
     *         {@code false} otherwise
     */
    public boolean isZeroMatrix() {
        return areAllEqual(0.0);
    }

    /**
     * Returns {@code true} if all entries of this {@code Matrix}
     * are approximately equal to the same constant value.
     *
     * <p>The first entry is used as the reference value. All other
     * entries are compared against it using the class's configured
     * numerical tolerance.
     *
     * <p><strong>Time Complexity:</strong> O(mn)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if all entries are approximately equal;
     *         {@code false} otherwise
     */
    public boolean isConstantMatrix() {
        return areAllEqual(entries[0][0]);
    }

    /**
     * Returns {@code true} if this {@code Matrix} is singular.
     *
     * <p>A matrix is singular if it is square and its determinant
     * is approximately equal to 0.0.
     *
     * <p>Determinant comparison uses the class's configured numerical
     * tolerance to account for floating-point rounding behavior.
     *
     * <p>If this matrix is not square, it is considered singular.
     *
     * <p><strong>Time Complexity:</strong> O(n³)<br>
     * <strong>Space Complexity:</strong> O(n²)
     *
     * @return {@code true} if this matrix is singular;
     *         {@code false} otherwise
     */
    public boolean isSingular() {
        return almostEqual(determinant(), 0.0);
    }

    /**
     * Returns {@code true} if this {@code Matrix} is symmetric.
     *
     * <p>A matrix is considered symmetric if:
     * <ul>
     *   <li>It is square.</li>
     *   <li>Each entry satisfies A[i][j] ≈ A[j][i].</li>
     * </ul>
     *
     * <p>Comparisons are performed using the class's configured
     * numerical tolerance to account for floating-point rounding behavior.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if this matrix is approximately symmetric;
     *         {@code false} otherwise
     */
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

    /**
     * Returns {@code true} if this {@code Matrix} is skew-symmetric.
     *
     * <p>A matrix is considered skew-symmetric if:
     * <ul>
     *   <li>It is square.</li>
     *   <li>Each entry satisfies A[i][j] ≈ -A[j][i].</li>
     *   <li>All diagonal entries are approximately equal to 0.0.</li>
     * </ul>
     *
     * <p>Comparisons are performed using the class's configured
     * numerical tolerance to account for floating-point rounding behavior.
     *
     * <p><strong>Time Complexity:</strong> O(n²)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return {@code true} if this matrix is approximately skew-symmetric;
     *         {@code false} otherwise
     */
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

    // ==== OBJECT METHODS ====

    /**
     * Returns a string representation of this {@code Matrix}.
     *
     * <p>The returned string presents the matrix in row-major order,
     * with each row on a separate line. Entries are formatted in a
     * readable, aligned manner to aid debugging and inspection.
     *
     * <p>This method is intended for diagnostic purposes and does
     * not guarantee a specific serialization format.
     *
     * <p><strong>Time Complexity:</strong> O(mn)<br>
     * <strong>Space Complexity:</strong> O(mn)
     *
     * @return a human-readable string representation of this matrix
     */
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

    /**
     * Returns a hash code for this {@code Matrix}.
     *
     * <p>The hash code is computed using the matrix dimensions
     * and the exact values of its entries (via
     * {@code Arrays.deepHashCode}).
     *
     * <p>This implementation is consistent with {@link #equals(Object)},
     * which performs strict element-wise comparison without tolerance.
     *
     * <p><strong>Time Complexity:</strong> O(mn)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @return a hash code value for this matrix
     */
    @Override
    public int hashCode() {
        return 31 * order.hashCode() + Arrays.deepHashCode(this.entries);
    }

    /**
     * Compares this {@code Matrix} with the specified object for equality.
     *
     * <p>Two matrices are considered equal if and only if:
     * <ul>
     *   <li>The specified object is also a {@code Matrix}.</li>
     *   <li>They have identical dimensions.</li>
     *   <li>All corresponding entries are exactly equal.</li>
     * </ul>
     *
     * <p>Entry comparison uses strict element-wise equality
     * (via {@code Double.compare}), not tolerance-based comparison.
     * Numerical tolerance is intentionally <strong>not</strong> applied
     * in order to preserve the general {@code equals}/{@code hashCode}
     * contract.
     *
     * <p><strong>Time Complexity:</strong> O(mn)<br>
     * <strong>Space Complexity:</strong> O(1)
     *
     * @param obj the object to compare with
     * @return {@code true} if the specified object represents a matrix
     *         with identical dimensions and entries; {@code false} otherwise
     */
    @Override
    public boolean equals(Object other) {
        if (!(other instanceof Matrix)) return false;

        Matrix o = (Matrix)other;
        return this.equalsMatrix(o);
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

    /**
     * Compares two double values using the configured equality tolerance.
     *
     * <p>This method is used for structural and numerical property
     * checks but is intentionally NOT used in {@link #equals(Object)}.
     *
     * @param a first value
     * @param b second value
     * @return {@code true} if the absolute difference between
     *         {@code a} and {@code b} is within the configured tolerance
     */
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

    /**
     * Returns a deep copy of the specified 2D array.
     *
     * <p>Each row is copied individually to ensure no aliasing
     * with the original grid.
     *
     * @param source the grid to copy
     * @return a deep copy of the provided grid
     */
    private double[][] deepGridCopy(double[][] grid) {

        int nrows = grid.length;
        int ncols = grid[0].length;
        double[][] result = new double[nrows][ncols];

        for (int i = 0; i < nrows; i++) {
            System.arraycopy(grid[i], 0, result[i], 0, ncols);
        }
        return result;
    }

    /**
     * Attempts to locate a valid pivot for the specified column
     * during row-reduction.
     *
     * <p>A pivot is considered valid if its absolute value exceeds
     * the configured pivot tolerance. This prevents division by
     * numerically unstable near-zero values.
     *
     * <p>If necessary, a row below the current pivot position may
     * be selected, indicating that a row swap is required.
     *
     * @param grid the working matrix during elimination
     * @param column the pivot column index
     * @return a {@code Pivot} describing the row index of the pivot
     *         and whether a row swap is required; returns a pivot
     *         with row = -1 if no valid pivot is found
     */
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

    /**
     * Swaps two rows in the provided grid in-place.
     *
     * <p>This method is used internally during elimination and
     * determinant computation to maintain correct pivot positioning.
     *
     * @param grid the matrix grid
     * @param r1 the first row index
     * @param r2 the second row index
     */
    private void swapGridRow(double[][] grid, int c, int row) {
        double[] temp = grid[c];
        grid[c] = grid[row];
        grid[row] = temp;
    }

    /**
     * Computes the product of the diagonal entries of the
     * specified square grid.
     *
     * <p>This method assumes the grid is square and does not
     * perform validation.
     *
     * @param grid a square matrix grid
     * @return the product of its diagonal elements
     */
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


}

