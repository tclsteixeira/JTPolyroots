using System;
using JTPolyroots;


namespace JTTests
{


    /// <summary>
    /// Implements tests for polynomials with real or complex (double precision) coeficents.
    /// </summary>
    /// <remarks>
    /// For polynomials with only real coeficients, use RPoly instead because is 4 times faster.
    /// </remarks>
    public static class TestCPoly
    {

        private const int C_FOR_START_IDX = 1;  // In Fortran array indexes are not zero based (start at 1)


        /// <summary>
        /// DRIVER TO TEST CPOLY.
        /// Tests the JT_CPoly Roots(...) procedure using some examples with known results.
        /// </summary>
        public static void Test()
        {
            DateTime starttime;
            //DateTime endtime;
            TimeSpan elapsedtime;
            JT_CPoly cpoly = new JT_CPoly();
            try
            {
                double[] _P = new double[14];
                double[] _PI = new double[14];
                double[] _ZR = new double[13];
                double[] _ZI = new double[13];
                bool _FAIL = false;
                string MSGFAIL = " CPOLY HAS FAILED ON THIS EXAMPLE";

                System.Console.WriteLine();
                System.Console.WriteLine("JENKINS AND TRAUB CPOLY TESTS DRIVER ****");
                System.Console.WriteLine("APPLIES TO POLYNOMIALS WITH COMPLEX AND REAL COEFICIENTS ****");
                System.Console.WriteLine();
                System.Console.WriteLine("EXAMPLE 1.  POLYNOMIAL WITH ZEROS 1,2,...,10.");
                System.Console.WriteLine();
                _P[1] = 1;
                _P[2] = -55;
                _P[3] = 1320;
                _P[4] = -18150;
                _P[5] = 157773;
                _P[6] = -902055;
                _P[7] = 3416930;
                _P[8] = -8409500;
                _P[9] = 12753576;
                _P[10] = -10628640;
                _P[11] = 3628800;

                for (int I = 1; I <= 11; I++)   //DO 10 I=1,11
                {
                    _PI[I] = 0.0;
                }

                PRTC(11, _P, _PI);

                starttime = DateTime.Now;

                // Solve poly
                cpoly.Roots(_P, _PI, 10, out _ZR, out _ZI, out _FAIL);

                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (_FAIL) goto _95;

                PRTZ(10, _ZR, _ZI);
                System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");
                System.Console.WriteLine();

            //----------------------------------------------
            _2: System.Console.WriteLine("EXAMPLE 2. ZEROS ON IMAGINARY AXIS DEGREE 3.");
                System.Console.WriteLine();

                _P[1] = 1;
                _P[2] = 0;
                _P[3] = -10001.0001;
                _P[4] = 0;
                _PI[1] = 0;
                _PI[2] = -10001.0001;
                _PI[3] = 0;
                _PI[4] = 1;
                PRTC(4, _P, _PI);

                starttime = DateTime.Now;

                // Solve poly
                cpoly.Roots(_P, _PI, 3, out _ZR, out _ZI, out _FAIL);

                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (_FAIL) goto _96;

                PRTZ(3, _ZR, _ZI);
                System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");
                System.Console.WriteLine();

            //----------------------------------------------
            _3: System.Console.WriteLine("EXAMPLE 3. ZEROS AT 1+I,1/2*(1+I)....1/(2**-9)*(1+I)");
                System.Console.WriteLine();

                _P[1] = 1.0;
                _P[2] = -1.998046875;
                _P[3] = 0.0;
                _P[4] = .7567065954208374;
                _P[5] = -.2002119533717632;
                _P[6] = 1.271507365163416E-2;
                _P[7] = 0;
                _P[8] = -1.154642632172909E-5;
                _P[9] = 1.584803612786345E-7;
                _P[10] = -4.652065399568528E-10;
                _P[11] = 0;
                _PI[1] = 0;
                _PI[2] = _P[2];
                _PI[3] = 2.658859252929688;
                _PI[4] = -7.567065954208374E-1;
                _PI[5] = 0;
                _PI[6] = _P[6];
                _PI[7] = -7.820779428584501E-4;
                _PI[8] = -_P[8];
                _PI[9] = 0;
                _PI[10] = _P[10];
                _PI[11] = 9.094947017729282E-13;
                PRTC(11, _P, _PI);

                starttime = DateTime.Now;

                // Solve poly
                cpoly.Roots(_P, _PI, 10, out _ZR, out _ZI, out _FAIL);

                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (_FAIL) goto _97;
                PRTZ(10, _ZR, _ZI);
                System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");
                System.Console.WriteLine();

            //----------------------------------------------
            _4: System.Console.WriteLine("EXAMPLE 4. MULTIPLE ZEROS");
                System.Console.WriteLine();
                _P[1] = 1;
                _P[2] = -10;
                _P[3] = 3;
                _P[4] = 284;
                _P[5] = -1293;
                _P[6] = 2374;
                _P[7] = -1587;
                _P[8] = -920;
                _P[9] = 2204;
                _P[10] = -1344;
                _P[11] = 288;
                _PI[1] = 0;
                _PI[2] = -10;
                _PI[3] = 100;
                _PI[4] = -334;
                _PI[5] = 200;
                _PI[6] = 1394;
                _PI[7] = -3836;
                _PI[8] = 4334;
                _PI[9] = -2352;
                _PI[10] = 504;
                _PI[11] = 0;
                PRTC(11, _P, _PI);

                starttime = DateTime.Now;

                // Solve poly
                cpoly.Roots(_P, _PI, 10, out _ZR, out _ZI, out _FAIL);

                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (_FAIL) goto _98;
                PRTZ(10, _ZR, _ZI);
                System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");
                System.Console.WriteLine();

            //----------------------------------------------
            _5: System.Console.WriteLine("EXAMPLE 5. 12 ZEROS EVENLY DISTRIBUTE ON A CIRCLE OF RADIUS 1. CENTERED AT 0+2I.");
                System.Console.WriteLine();
                _P[1] = 1;
                _P[2] = 0;
                _P[3] = -264;
                _P[4] = 0;
                _P[5] = 7920;
                _P[6] = 0;
                _P[7] = -59136;
                _P[8] = 0;
                _P[9] = 126720;
                _P[10] = 0;
                _P[11] = -67584;
                _P[12] = 0;
                _P[13] = 4095;
                _PI[1] = 0;
                _PI[2] = -24;
                _PI[3] = 0;
                _PI[4] = 1760;
                _PI[5] = 0;
                _PI[6] = -25344;
                _PI[7] = 0;
                _PI[8] = 101376;
                _PI[9] = 0;
                _PI[10] = -112640;
                _PI[11] = 0;
                _PI[12] = 24576;
                _PI[13] = 0;
                PRTC(13, _P, _PI);

                starttime = DateTime.Now;
                // Solve poly
                cpoly.Roots(_P, _PI, 12, out _ZR, out _ZI, out _FAIL);
                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (_FAIL) goto _99;

                PRTZ(12, _ZR, _ZI);
                System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");

                System.Console.WriteLine();
                System.Console.WriteLine("EXAMPLE 6. POLYNOMIAL WITH ZEROS 1,2,...,10.'");
                System.Console.WriteLine();

                int degree = 10;
                _P = new double[degree + 2];
                _PI = new double[degree + 2];
                _PI.Initialize();

                _P[1] = 1.0;
                _P[2] = -55.0;
                _P[3] = 1320.0;
                _P[4] = -18150.0;
                _P[5] = 157773.0;
                _P[6] = -902055.0;
                _P[7] = 3416930.0;
                _P[8] = -8409500.0;
                _P[9] = 12753576.0;
                _P[10] = -10628640.0;
                _P[11] = 3628800.0;

                //Print coeficients.
                PRTC(degree + 1, _P, _PI);

                starttime = DateTime.Now;

                //Solve poly.
                cpoly.Roots(_P, _PI, degree, out _ZR, out _ZI, out _FAIL);

                elapsedtime = TimeSpan.FromTicks(DateTime.Now.Ticks - starttime.Ticks);

                if (_FAIL)
                {
                    System.Console.WriteLine(" ** Failure by CPOLY **");
                }
                else
                {
                    //Output zeros.
                    PRTZ(degree, _ZR, _ZI);
                    System.Console.WriteLine("Elapsed time: " + elapsedtime.Milliseconds + " ms");
                }

                System.Console.WriteLine();
                System.Console.WriteLine("END CPOLY TESTS ***************");
                System.Console.WriteLine();
                System.Console.WriteLine("press ENTER to continue with RPOLY tests...");
                System.Console.ReadLine();
                return;

            _95: System.Console.WriteLine(MSGFAIL);
                goto _2;
            _96: System.Console.WriteLine(MSGFAIL);
                goto _3;
            _97: System.Console.WriteLine(MSGFAIL);
                goto _4;
            _98: System.Console.WriteLine(MSGFAIL);
                goto _5;
            _99: System.Console.WriteLine(MSGFAIL);
                System.Console.ReadLine();
                return;
            }
            catch (Exception ex)
            {
                System.Console.WriteLine(ex.ToString());
                System.Console.ReadLine();
                throw ex;
            }
        }


        /// <summary>
        /// Prints the polinomial coeficients.
        /// </summary>
        /// <param name="N">Number of coeficients.</param>
        /// <param name="P">The real part coeficients array.</param>
        /// <param name="Q">The imaginary part coeficients array.</param>
        private static void PRTC(int N, double[] P, double[] Q)
        {
            try
            {
                for (int I = C_FOR_START_IDX; I <= N; I++)
                {
                    System.Console.WriteLine(string.Format("COEFICIENTS: {0}, {1}", P[I], Q[I]));
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
                    System.Console.WriteLine(string.Format("ZEROS: {0} {1}", ZR[I], ZI[I]));
                }
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }



    }


}
