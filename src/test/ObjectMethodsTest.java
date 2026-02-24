package test;

import main.Matrix;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

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
                new double[] { 0, -7 },
                new double[] { 123456, -999999 }
        );

        A2 = Matrix.ofRows(
                new double[] { 0, -7 },
                new double[] { 123456, -999999 }
        );

        B1 = Matrix.ofRows(
                new double[] { 1.0000001, -2000.5, 0 },
                new double[] { 3.1415926, 999999.999, -42.42 },
                new double[] { -0.0000003, 7, -123456.789 }
        );

        B2 = Matrix.ofRows(
                new double[] { 1.0000002, -2000.5000004, 0.0000001 },
                new double[] { 3.1415925, 999999.9990003, -42.4200002 },
                new double[] { -0.0000004, 7.0000003, -123456.7889996 }
        );

        C1 = Matrix.ofRows(
                new double[] { 1, -2, 3, -4 },
                new double[] { 1000.5, 0, -999.25, 42 },
                new double[] { -7, 8.8, 0.0001, -50000 }
        );

        C2 = Matrix.ofRows(
                new double[] { 2, -3, 4, -5 },
                new double[] { 1001.5, 1, -998.25, 43 },
                new double[] { -6, 9.8, 1.0001, -49999 }
        );

        U1 = Matrix.ofRows(
                new double[] { 5.0 }
        );

       U2 = Matrix.ofRows(
                new double[] { 5.0 + 0.9e-6 }
        );

       D1 = Matrix.ofRows(
                new double[] { 10.0 }
        );

       D2 = Matrix.ofRows(
                new double[] { 10.0 + 1e-6 }
        );

       O1 = Matrix.ofRows(
                new double[] { -3.0 }
        );

       O2 = Matrix.ofRows(
                new double[] { -3.0 + 1.1e-6 }
        );

    }


    @Test
    @DisplayName("Equals method works.")
    void doNearlyEqualMatricesReturnTrue() {

        //reflexivity
        assertEquals(A1, A1);
        assertEquals(A2, A2);
        assertEquals(B1, B1);
        assertEquals(B2, B2);
        assertEquals(C1, C1);
        assertEquals(C2, C2);

        //symmetry
        assertEquals(A1, A2);
        assertEquals(A2, A1);
        assertNotEquals(B1, B2);
        assertNotEquals(B2, B1);

        //null checks
        assertNotEquals(null, A1);
        assertNotEquals("This is a Matrix", A2);

        //dimension-based inequality
        assertNotEquals(new Matrix(3,2), new Matrix(2,3));
        assertNotEquals(new Matrix(5,5), new Matrix(5,3));

        //standard comparisons
        assertEquals(A1, A2);
        assertNotEquals(A1, B1);
        assertNotEquals(A1, B2);

        assertNotEquals(B1, B2);
        assertNotEquals(B1, C1);
        assertNotEquals(B1, C2);
        assertNotEquals(B2, C1);
        assertNotEquals(B2, C2);

        assertNotEquals(C1, C2);

        //Edge cases
        assertNotEquals(U1, U2);
        assertNotEquals(D1, D2);
        assertNotEquals(O1, O2);
        
    }

    @Test
    @DisplayName("Hashcode method works.")
    void doNearlyEqualMatricesHaveSimilarHashcodes() {

        assertEquals(A1.hashCode(), A2.hashCode());
        assertNotEquals(A1.hashCode(), B1.hashCode());
        assertNotEquals(A1.hashCode(), B2.hashCode());

        assertNotEquals(B1.hashCode(), B2.hashCode());
        assertNotEquals(B1.hashCode(), C1.hashCode());
        assertNotEquals(B1.hashCode(), C2.hashCode());
        assertNotEquals(B2.hashCode(), C1.hashCode());
        assertNotEquals(B2.hashCode(), C2.hashCode());

        assertNotEquals(C1.hashCode(), C2.hashCode());
        assertNotEquals(U1.hashCode(), U2.hashCode());
        assertNotEquals(D1.hashCode(), D2.hashCode());
        assertNotEquals(D1.hashCode(), D2.hashCode());
    }

    @Test
    @DisplayName("Class conforms to equals-hashcode contract.")
    void doesMatrixConformToEqualsHashcodeContract() {

        if (A1.equals(A2)) {
            assertEquals(A1.hashCode(), A2.hashCode());
        }

        if (B1.equals(B2)) {
            assertEquals(B1.hashCode(), B2.hashCode());
        }

    }

    @Test
    void doesToStringHandleEntriesWithDiverseSizes() {

    }
}
