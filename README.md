# Java Matrix Library by Arnav Merani

A numerically aware, row-major Java implementation of a fully functional matrix data structure supporting core linear algebra operations, structural queries, and explicit equality semantics. A high-performance, Maven-managed Linear Algebra library for Java 21 including end-to-end Javadoc and wide JUnit testing coverage.

---

## Features

- Matrix addition and subtraction
- Matrix multiplication
- Scalar multiplication
- Determinant computation (Gaussian elimination-based)
- Matrix inverse (Gauss–Jordan elimination with partial pivoting)
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

### 1. Install the Library Locally

Clone the repository and install the artifact into your local Maven repository:

```bash
   mvn clean install
```
This installs the library into:

```bash
  ~/.m2/repository/com/arnavmerani/matrix/1.0.1/
```

### 2. Add dependency to your project

Create a new project with the Maven build system, and add the
following dependency to your pom.xml file:

```XML
<dependency>
    <groupId>com.arnavmerani</groupId>
    <artifactId>matrix</artifactId>
    <version>1.0.1</version>
</dependency>
```

Then reload Maven in your project using:
```bash
  mvn clean compile 
```


## Usage Example
```Java
import com.arnavmerani.matrix.Matrix;

public class Main {
    public static void main(String[] args) {

        Matrix A = new Matrix(new double[][] {
            {1, 2},
            {3, 4}
        });

        Matrix B = A.inverse();
        double det = A.determinant();

        Matrix C = A.multiply(B);

        System.out.println("det(A) = " + det);
        System.out.println("A * A⁻¹ =");
        System.out.println(C);
    }
}
```
---

## License and Contact
- Maintainer: Arnav Merani
- GitHub: arnavmer-935
- License: MIT
