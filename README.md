# Matrix.java

A Java implementation of a fully functional matrix data structure supporting core linear algebra operations, including arithmetic, determinant computation, matrix inversion, and Gauss–Jordan row reduction.

---

## Features

- Matrix addition and subtraction  
- Matrix multiplication  
- Scalar multiplication  
- Determinant computation  
- Matrix inverse (via Gauss–Jordan elimination)  
- Row reduction to Reduced Row Echelon Form (RREF)  
- In-place and out-of-place operation variants  
- Strict dimension and argument validation  
- Floating-point tolerance handling (~1e-6) for numerical stability  

---

## Implementation Highlights

- Determinant and inverse implemented using systematic row operations  
- Pivot selection and row-swapping logic for stable elimination  
- Defensive programming for singular and non-square matrices  
- Clear API design separating mutating and non-mutating methods  
- Emphasis on correctness and edge-case handling  

---

## Testing

- Validation against singular matrices and invalid dimensions  
- Algebraic identity checks where applicable  
- Incremental expansion of test coverage  
- Numerical comparisons performed within defined tolerance bounds  

---

## Technical Focus

This project emphasizes:

- Algorithm implementation from first principles  
- Numerical reasoning and floating-point precision management  
- Object-oriented design and API structuring  
- Robust input validation and error handling  

---

## Motivation

Developed to deepen understanding of linear algebra concepts through direct implementation, translating mathematical definitions into reliable, well-structured code.
