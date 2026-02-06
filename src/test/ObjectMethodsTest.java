package test;

import main.Matrix;
import main.MatrixException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static java.lang.Double.NaN;
import static org.junit.jupiter.api.Assertions.*;

public class ObjectMethodsTest {
    
    private static final double TOLERANCE = 1e-6;
    private Matrix A1, A2;
    private Matrix B1, B2;
    private Matrix C1, C2;

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
    }


    @Test
    void doNearlyEqualMatricesReturnTrue() {

        assertEquals(A1, A2);
        assertNotEquals(A1, B1);
        assertNotEquals(A1, B2);

        assertEquals(B1, B2);
        assertNotEquals(B1, C1);
        assertNotEquals(B1, C2);
        assertNotEquals(B2, C1);
        assertNotEquals(B2, C2);

        assertNotEquals(C1, C2);
        
    }

    @Test
    void doNearlyEqualMatricesHaveSimilarHashcodes() {

        assertEquals(A1.hashCode(), A2.hashCode(), TOLERANCE);
        assertNotEquals(A1.hashCode(), B1.hashCode(), TOLERANCE);
        assertNotEquals(A1.hashCode(), B2.hashCode(), TOLERANCE);

        assertEquals(B1.hashCode(), B2.hashCode(), TOLERANCE);
        assertNotEquals(B1.hashCode(), C1.hashCode(), TOLERANCE);
        assertNotEquals(B1.hashCode(), C2.hashCode(), TOLERANCE);
        assertNotEquals(B2.hashCode(), C1.hashCode(), TOLERANCE);
        assertNotEquals(B2.hashCode(), C2.hashCode(), TOLERANCE);

        assertNotEquals(C1.hashCode(), C2.hashCode(), TOLERANCE);
    }

    @Test
    void doesToStringHandleEntriesWithDiverseSizes() {

    }
}
