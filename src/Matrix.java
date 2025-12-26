import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class Matrix {

    // ==== INSTANCE VARIABLES AND DECLARATIONS ====
    private record Pair(int x, int y) {
        @Override
        public String toString() {
            return String.format("(%d, %d)", x, y);
        }
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

        if (grid.length == 0 || grid[0].length == 0) {
            throw MatrixException.illegalDimensions();
        }

        if (isJaggedGrid(grid)) {
            throw MatrixException.jaggedMatrix(getMismatchedRowIndex(grid));
        }

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
            for (int j = 0; j < columns; j++) {
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
        if (rows == null) {
            throw new IllegalArgumentException("Matrix grid must be non-null.");
        }

        if (rows.length <= 0 || rows[0].length == 0) {
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

    public static Matrix nullMatrix(int nrows, int ncols) {
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

    public void setEntry(double value, int r, int c) {
        if (!isInBounds(r,c)) {
            throw new IndexOutOfBoundsException(String.format("Cell coordinates (%d, %d) out of bounds for Matrix of order %s", r, c, order));
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

        this.entries[rowIndex] = row.clone();
    }

    public void setColumn(double[] col, int colIndex) {
        if (col.length != rows) {
            throw MatrixException.columnLengthMismatch();

        }

        if (!colInRange(colIndex)) {
            throw new IndexOutOfBoundsException(String.format("Column index %d is out of bounds for Matrix of order %s", colIndex, this.order));
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
            return nullMatrix(this.rows, this.columns);
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
            return nullMatrix(this.rows, other.columns);
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

    public Matrix inverse(Matrix m) {
        if (isSingular(m)) {
            throw MatrixException.matrixSingularity();
        }

        int nrows = m.getRows();
        int ncols = m.getColumns();
        double[][] adjoint = new double [nrows][ncols];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; ++j) {
                adjoint[j][i] = determinant(getCofactorMatrix(i, j));
            }
        }

        Matrix adj = new Matrix(adjoint);
        return adj.multiplyByScalar(Math.pow(determinant(m), -1));
    }

    // ==== NUMERICAL METHODS ====
    //TODO: Optimize using row reduction
    public double determinant(Matrix m) {
        if (!m.isSquareMatrix()) {
            return 0;
        }

        if (m.getOrder().equals(new Pair(1,1))) { //single element matrix
            return this.entries[0][0];
        }

        if (m.getOrder().equals(new Pair(2,2))) { //2x2 square matrix
            return (m.getValue(0,0) * m.getValue(1,1)) - (m.getValue(0,1) * m.getValue(1,0));
        }

        double[] topRow = m.entries[0];
        double det = 0;
        for (int j = 0; j < topRow.length; j++) {
            det += ((int)Math.pow(-1, j) * getValue(0,j) * determinant(m.getCofactorMatrix(0,j)));
            //determinant of a matrix is sum of products of elements in a row/column times that element's cofactor
            //cofactor = -1^(i+j) * determinant of the minor of
        }
        return det;
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

    public boolean isSingular(Matrix m) {
        return almostEqual(determinant(m), 0.0);
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
        return 0 <= r && r < this.rows && 0 <= c && c < this.columns;
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

    private Matrix getCofactorMatrix(int mthRow, int nthCol) {
        if (!isInBounds(mthRow, nthCol)) {
            throw new IndexOutOfBoundsException(String.format("Cell coordinates (%d, %d) out of bounds for Matrix of order %s", mthRow, nthCol, order));
        }

        List<List<Double>> gridList = new ArrayList<>(); //convert to arraylist grid for easier deletion of row and column
        for (double[] row : this.entries) {
            gridList.add(Arrays.stream(row).boxed().collect(Collectors.toList()));
        }

        gridList.forEach(rowList -> {
            rowList.remove(nthCol);
        });

        gridList.remove(mthRow);

        return new Matrix(gridList);
    }

    private double[][] populateJaggedMatrix(double[][] matrix) {
        if (matrix.length == 0) {
            throw new IllegalArgumentException("Matrix grid must be non-empty");
        }

        double[] longestRow = {};
        for (double[] row : matrix) {
            if (longestRow.length < row.length) {
                longestRow = row;
            }
        }

        double[][] populatedMatrix = new double[matrix.length][longestRow.length];
        for (int i = 0; i < matrix.length; i++) {
            System.arraycopy(matrix[i], 0, populatedMatrix[i], 0, matrix[i].length);
        }

        return populatedMatrix;
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
            for (int j = 0; j < ncols; j++) {
                result[i][j] = grid[i][j];
            }
        }
        return result;
    }

    private String formattedValue(double v) {
        return Math.abs(v) <= TOLERANCE ? "0.0" : String.format("%.3f", v);
    }

    // ==== OBJECT METHODS ====

    @Override
    public String toString() {
        StringBuilder matrix = new StringBuilder();
        for (int r = 0; r < rows; r++) {
            matrix.append("[ ");
            for (int c = 0; c < columns; c++) {
                matrix.append(formattedValue(getValue(r,c)));
                matrix.append(" ");
            }
            matrix.append("]\n");
        }
        return matrix.toString();
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

