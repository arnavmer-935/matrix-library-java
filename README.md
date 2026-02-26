# Matrix.java

A numerically aware, row-major Java implementation of a fully functional matrix data structure supporting core linear algebra operations, structural queries, and explicit equality semantics.

---

## Features

- Matrix addition and subtraction
- Matrix multiplication
- Scalar multiplication
- Determinant computation (Gaussian elimination-based)
- Matrix inverse (Gauss–Jordan elimination with partial pivoting)
- Row reduction to Reduced Row Echelon Form (RREF)
- Symmetric and skew-symmetric decomposition
- Structural queries (square, triangular, symmetric, identity, etc.)
- In-place and out-of-place operation variants
- Strict dimension and argument validation

---

## Numerical Design

- Structural and property checks use a configured floating-point tolerance (~1e-6) to account for rounding behavior.
- Pivot detection during elimination uses a stricter tolerance to handle near-singular matrices.
- `equals()` performs strict element-wise comparison using `Double.compare` (no tolerance).
- `hashCode()` is consistent with strict equality semantics.

This separation ensures numerical robustness while preserving the `equals`/`hashCode` contract.

---

## Mutability Model

- Methods ending with `InPlace` mutate the current instance.
- Non in-place methods return new `Matrix` instances and leave the original unchanged.
- Internal storage is deeply owned; external mutation is prevented via defensive copying.

---

## Getting Started

```java

import com.arnavmerani.matrix

Matrix A = Matrix.ofRows(new double[][] {
    {1, 2},
    {3, 4}
});

Matrix B = A.inverse();
double det = A.determinant();

Matrix C = A.multiply(B);
System.out.println(C);