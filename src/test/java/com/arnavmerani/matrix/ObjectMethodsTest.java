package com.arnavmerani.matrix;

import com.arnavmerani.matrix.Matrix;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class ObjectMethodsTest {

    private Matrix A1, A2;
    private Matrix B1, B2;
    private Matrix C1, C2;
    private Matrix U1, U2;
    private Matrix D1, D2;
    private Matrix O1, O2;

    @BeforeEach
    void setup() {
        A1 = Matrix.ofRows(
                new double[]{0, -7},
                new double[]{123456, -999999}
        );

        A2 = Matrix.ofRows(
                new double[]{0, -7},
                new double[]{123456, -999999}
        );

        B1 = Matrix.ofRows(
                new double[]{1.0000001, -2000.5, 0},
                new double[]{3.1415926, 999999.999, -42.42},
                new double[]{-0.0000003, 7, -123456.789}
        );

        B2 = Matrix.ofRows(
                new double[]{1.0000002, -2000.5000004, 0.0000001},
                new double[]{3.1415925, 999999.9990003, -42.4200002},
                new double[]{-0.0000004, 7.0000003, -123456.7889996}
        );

        C1 = Matrix.ofRows(
                new double[]{1, -2, 3, -4},
                new double[]{1000.5, 0, -999.25, 42},
                new double[]{-7, 8.8, 0.0001, -50000}
        );

        C2 = Matrix.ofRows(
                new double[]{2, -3, 4, -5},
                new double[]{1001.5, 1, -998.25, 43},
                new double[]{-6, 9.8, 1.0001, -49999}
        );

        U1 = Matrix.ofRows(new double[]{5.0});
        U2 = Matrix.ofRows(new double[]{5.0 + 0.9e-6});

        D1 = Matrix.ofRows(new double[]{10.0});
        D2 = Matrix.ofRows(new double[]{10.0 + 1e-6});

        O1 = Matrix.ofRows(new double[]{-3.0});
        O2 = Matrix.ofRows(new double[]{-3.0 + 1.1e-6});
    }

    @Nested
    class EqualsTests {

        @Test
        @DisplayName("Reflexivity: matrix equals itself")
        void reflexivity() {
            assertEquals(A1, A1);
            assertEquals(B1, B1);
            assertEquals(C1, C1);
        }

        @Test
        @DisplayName("Symmetry for equal and unequal matrices")
        void symmetry() {
            assertEquals(A1, A2);
            assertEquals(A2, A1);

            assertNotEquals(B1, B2);
            assertNotEquals(B2, B1);
        }

        @Test
        @DisplayName("Handles null and different object types")
        void nullAndTypeChecks() {
            assertNotEquals(null, A1);
            assertNotEquals("not a matrix", A1);
        }

        @Test
        @DisplayName("Dimension mismatch implies inequality")
        void dimensionMismatch() {
            assertNotEquals(new Matrix(3,2), new Matrix(2,3));
            assertNotEquals(new Matrix(5,5), new Matrix(5,3));
        }

        @Test
        @DisplayName("Distinct matrices are not equal")
        void distinctMatrices() {
            assertNotEquals(A1, B1);
            assertNotEquals(B1, C1);
            assertNotEquals(C1, C2);
        }

        @Test
        @DisplayName("Strict equality: near-equal doubles are not equal")
        void strictDoubleEquality() {
            assertNotEquals(U1, U2);
            assertNotEquals(D1, D2);
            assertNotEquals(O1, O2);
        }
    }

    @Nested
    class HashCodeTests {

        @Test
        @DisplayName("Equal matrices have equal hash codes")
        void equalMatrices_sameHash() {
            assertEquals(A1.hashCode(), A2.hashCode());
        }

        @Test
        @DisplayName("Unequal matrices typically have different hash codes")
        void unequalMatrices_differentHash() {
            assertNotEquals(A1.hashCode(), B1.hashCode());
            assertNotEquals(B1.hashCode(), C1.hashCode());
            assertNotEquals(C1.hashCode(), C2.hashCode());
            assertNotEquals(U1.hashCode(), U2.hashCode());
        }
    }

    @Nested
    class EqualsHashCodeContract {

        @Test
        @DisplayName("If equals returns true, hash codes must match")
        void equalsImpliesSameHash() {
            if (A1.equals(A2)) {
                assertEquals(A1.hashCode(), A2.hashCode());
            }
        }
    }
}