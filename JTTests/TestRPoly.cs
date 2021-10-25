using System;
using JTPolyroots;

namespace JTTests
{

    /// <summary>
    /// Implements tests for polynomials with real (double precision) coeficents.
    /// </summary>
    public static class TestRPoly
    {

        private const int C_FOR_START_IDX = 1;  // In Fortran array indexes are not zero based (start at 1)


        /// <summary>
        /// Tests the JT_RPoly Roots(...) procedure using some examples with known results.
        /// </summary>
        /// <remarks>
        /// JT_RPoly Roots(...) solves polynomials with real (double precision) coeficients
        /// to find its roots.
        /// </remarks>
        public static void Test()
        {
            DateTime starttime;
            TimeSpan elapsedtime;
            try
            {
                JT_RPoly rpoly = new JT_RPoly();

                double[] p, zr, zi;
                int degree;//, i;
                bool fail;

                System.Console.WriteLine();
                System.Console.WriteLine("JENKINS AND TRAUB RPOLY TESTS DRIVER ****");
                System.Console.WriteLine("APPLIES ONLY TO POLYNOMIALS WITH REAL COEFICIENTS ****");
                System.Console.WriteLine();
                WRITE("EXAMPLE 1. POLYNOMIAL WITH ZEROS 1,2,...,10.'");
                System.Console.WriteLine();

                degree = 10;
                p = new double[degree + 2];

                p[1] = 1.0;
                p[2] = -55.0;
                p[3] = 1320.0;
                p[4] = -18150.0;
                p[5] = 157773.0;
                p[6] = -902055.0;
                p[7] = 3416930.0;
                p[8] = -8409500.0;
                p[9] = 12753576.0;
                p[10] = -10628640.0;
                p[11] = 3628800.0;

                // Print coeficients.
                PRTC(11, p);

                starttime = DateTime.Now;

                // Solve poly.
                rpoly.Roots(p, degree, out zr, out zi, out fail);

                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (fail)
                {
                    WRITE(" ** Failure by RPOLY **");
                }
                else
                {
                    //Output zeros.
                    PRTZ(degree, zr, zi);
                    System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");
                    System.Console.WriteLine();
                }

                //! This test provided by Larry Wigton
                WRITE("");
                WRITE("Now try case where 1 is an obvious root");
                WRITE("");

                degree = 5;
                p = new double[degree + 2];

                p[1] = 8.0;
                p[2] = -8.0;
                p[3] = 16.0;
                p[4] = -16.0;
                p[5] = 8.0;
                p[6] = -8.0;

                //Print coeficients.
                PRTC(6, p);

                starttime = DateTime.Now;

                //Solve poly.
                rpoly.Roots(p, degree, out zr, out zi, out fail);

                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (fail)
                {
                    WRITE(" ** Failure by RPOLY **");
                }
                else
                {
                    PRTZ(degree, zr, zi);
                    System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");
                    System.Console.WriteLine();
                }

                WRITE("");
                System.Console.WriteLine("END RPOLY TESTS ***************");
                WRITE("");
                WRITE("Press ENTER to exit...");
                READLN();
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Writes to the output screen the specified text.
        /// </summary>
        /// <param name="str">The source text to write.</param>
        private static void WRITE(string txt)
        {
            try
            {
                Console.WriteLine(txt);
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Prints the polinomial coeficients.
        /// </summary>
        /// <param name="N">Number of coeficients.</param>
        /// <param name="P">The real part coeficients array.</param>
        private static void PRTC(int N, double[] P)    //, double[] Q) /// <param name="Q">The imaginary part coeficients array.</param>
        {
            try
            {
                for (int I = C_FOR_START_IDX; I <= N; I++)
                {
                    System.Console.WriteLine(string.Format("COEFICIENTS: {0}", P[I]));  //, Q[I]));
                }
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Prints the polynomial zeros (roots).
        /// </summary>
        /// <param name="N">Number of zeros.</param>
        /// <param name="ZR">Real part zeros array.</param>
        /// <param name="ZI">Imaginary part zeros array.</param>
        private static void PRTZ(int N, double[] ZR, double[] ZI)
        {
            try
            {
                for (int I = C_FOR_START_IDX; I <= N; I++)
                {
                    WRITE(string.Format("ZEROS: {0} {1}", ZR[I], ZI[I]));
                }
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Waits for keyboard carrige return press.
        /// </summary>
        private static void READLN()
        {
            try
            {
                Console.ReadLine();
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }

    }
}
