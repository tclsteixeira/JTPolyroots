using System;

namespace JTPolyroots
{

    /// <summary>
    /// Jenkins-Traub complex (and real) polynomial root finder.
    /// --- Applies only to one variable polynomials (with real and complex coeficients) ---
    /// </summary>
    /// <remarks>
    /// This is a direct port from FORTRAN77.
    /// It works for polynomials with Real and Complex coeficients, but if you are using only 
    /// polynomials with Real coeficients use instead JT_RPoly because it is 4 times faster.
    /// It's one of the most valuable method for polynomial root finding used for ex. in Maple and Mathematica soft.
    /// 
    ///      ALGORITHM 419 COLLECTED ALGORITHMS FROM ACM.
    ///      ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 02,
    ///      P. 097.
    /// 
    ///      Added changes from Remark on Algorithm 419 by David H. Withers
    ///      CACM (March 1974) Vol 17 No 3 p. 157 
    /// 
    /// </remarks>
    public sealed class JT_CPoly
    {

        private const int C_MAXDEEGRE = 50;
        private static int mMaxSize = C_MAXDEEGRE + 1;

        #region COMMON/GLOBAL/
        
        //Remmember that in FORTRAN array indexes start at 1.
        double[] PR = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] PI = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] HR = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] HI = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] QPR = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] QPI = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] QHR = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] QHI = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] SHR = new double[mMaxSize];//[C_MAXDEEGRE];
        double[] SHI = new double[mMaxSize];//[C_MAXDEEGRE];

        double SR, SI, TR, TI, PVR, PVI, ARE, MRE, ETA, INFIN;
        int NN;

        #endregion COMMON/GLOBAL/


        //Start index for fortran arrays.
        private const int C_FOR_START_IDX = 1;



        public JT_CPoly()
            : base()
        { }



        #region Private Methods


        /// <summary>
        ///     MCON PROVIDES MACHINE CONSTANTS USED IN VARIOUS PARTS OF THE
        ///     PROGRAM. THE USER MAY EITHER SET THEM DIRECTLY OR USE THE
        ///     STATEMENTS BELOW TO COMPUTE THEM. THE MEANING OF THE FOUR
        ///     CONSTANTS ARE -
        ///     ETA       THE MAXIMUM RELATIVE REPRESENTATION ERROR
        ///     WHICH CAN BE DESCRIBED AS THE SMALLEST POSITIVE
        ///     FLOATING-POINT NUMBER SUCH THAT 1.0D0 + ETA IS
        ///     GREATER THAN 1.0D0.
        ///     INFINY    THE LARGEST FLOATING-POINT NUMBER
        ///     SMALNO    THE SMALLEST POSITIVE FLOATING-POINT NUMBER
        ///     BASE      THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED
        ///     LET T BE THE NUMBER OF BASE-DIGITS IN EACH FLOATING-POINT
        ///     NUMBER(DOUBLE PRECISION). THEN ETA IS EITHER .5*B**(1-T)
        ///     OR B**(1-T) DEPENDING ON WHETHER ROUNDING OR TRUNCATION
        ///     IS USED.
        ///     LET M BE THE LARGEST EXPONENT AND N THE SMALLEST EXPONENT
        ///     IN THE NUMBER SYSTEM. THEN INFINY IS (1-BASE**(-T))*BASE**M
        ///     AND SMALNO IS BASE**N.
        ///     THE VALUES FOR BASE,T,M,N BELOW CORRESPOND TO THE IBM/360.
        /// </summary>
        /// <param name="ETA">
        ///     THE MAXIMUM RELATIVE REPRESENTATION ERROR
        ///     WHICH CAN BE DESCRIBED AS THE SMALLEST POSITIVE
        ///     FLOATING-POINT NUMBER SUCH THAT 1.0D0 + ETA IS
        ///     GREATER THAN 1.0D0.
        /// </param>
        /// <param name="INFINY">THE LARGEST FLOATING-POINT NUMBER</param>
        /// <param name="SMALNO">THE SMALLEST POSITIVE FLOATING-POINT NUMBER</param>
        /// <param name="_BASE">THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED</param>
        private void MCON(out double ETA, out double INFINY, out double SMALNO, out double _BASE)
        {
            try
            {
                int M, N, T;

                _BASE = 16.0E0;
                T = 14;
                M = 63;
                N = -65;
                ETA = Math.Pow(_BASE, 1 - T);
                INFINY = _BASE * (1.0E0 - Math.Pow(_BASE, -T) ) * Math.Pow(_BASE, M - 1);
                SMALNO = Math.Pow(_BASE, N + 3) / Math.Pow(_BASE, 3);
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        ///  COMPUTES THE MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.
        /// </summary>
        /// <param name="R">Complex number real part.</param>
        /// <param name="I">complex number imaginary part.</param>
        /// <returns>Returns the modulus of a complex number with a specified real and imaginary parts.</returns>
        private double CMOD(double R, double I)
        {
            double FResult = double.NaN;
            try 
	        {	        
		        double AR, AI;
                AR = Math.Abs(R);
                AI = Math.Abs(I);

                if (AR >= AI) goto _10;

                FResult = AI * Math.Sqrt(1.0 + Math.Pow(AR / AI , 2));
                return FResult;

          _10:  if (AR <= AI) goto _20;

                FResult = AR * Math.Sqrt(1.0 + Math.Pow(AI / AR, 2));
                return FResult;

          _20:  FResult = AR * Math.Sqrt( 2.0 );
                return FResult;
	        }
	        catch (Exception ex)
	        {
		        throw ex;
	        }

            //return FResult;
        }


        /// <summary>
        ///     RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE
        ///     POLYNOMIAL. THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID
        ///     UNDETECTED UNDERFLOW INTERFERING WITH THE CONVERGENCE
        ///     CRITERION.  THE FACTOR IS A POWER OF THE BASE.
        ///     PT - MODULUS OF COEFFICIENTS OF P
        ///     ETA,INFIN,SMALNO,BASE - CONSTANTS DESCRIBING THE
        ///     FLOATING POINT ARITHMETIC.
        /// </summary>
        /// <param name="NN"></param>
        /// <param name="PT">MODULUS OF COEFFICIENTS OF P.</param>
        /// <param name="ETA">CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC (epsylon).</param>
        /// <param name="INFIN">CONSTANTS DESCRIBING THE MAXIMUM FLOATING POINT ARITHMETIC VALUE.</param>
        /// <param name="SMALNO">CONSTANTS DESCRIBING THE MINIMUM FLOATING POINT ARITHMETIC VALUE.</param>
        /// <param name="_BASE">CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC BASE.</param>
        /// <returns>RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL</returns>
        private double SCALE(int NN, double[] PT, double ETA, double INFIN, double SMALNO, double _BASE)
        {
            double FResult = double.NaN;
            try
            {
                double HI, LO, MAX, MIN, X, SC, L;  //TODO: Check if L should be global.

                // FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
                HI = Math.Sqrt(INFIN);
                LO = SMALNO / ETA;
                MAX = 0.0;
                MIN = INFIN;

                for (int I = 1; I <= NN; I++)    //DO 10 I = 1,NN
			    {
                    X = PT[I];
                    if (X > MAX) MAX = X;
                    if ((X != 0.0) && (X < MIN)) MIN = X;
			    }

                // SCALE ONLY IF THERE ARE VERY LARGE OR VERY SMALL COMPONENTS.
                FResult = 1.0;
                if ((MIN >= LO) && (MAX <= HI)) return FResult;

                X = LO / MIN;
                if (X > 1.0) goto _20;

                SC = 1.0 / (Math.Sqrt(MAX) * Math.Sqrt(MIN));
                goto _30;

          _20:  SC = X;
                if ((INFIN / SC) > MAX) SC = 1.0;

          _30:  L = Math.Log(SC) / Math.Log(_BASE) + 0.500;
                FResult = Math.Pow(_BASE, L);
                return FResult;
            }
            catch (Exception ex)
            {
                throw ex;
            }

            //return FResult;
        }


        /// <summary>
        /// COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
        /// </summary>
        /// <param name="AR">Complex number A real part.</param>
        /// <param name="AI">Complex number A imaginary part.</param>
        /// <param name="BR">Complex number B real part.</param>
        /// <param name="BI">Complex number B imaginary part.</param>
        /// <param name="CR">Returns the complex number result real part.</param>
        /// <param name="CI">Returns the complex number result imaginary part.</param>
        private void CDIVID(double AR, double AI, double BR, double BI, out double CR, out double CI)
        {
            try
            {
                double R, D, T, INFIN;

                if ((BR != 0.0) || (BI != 0.0)) goto _10;

                // DIVISION BY ZERO, C = INFINITY.
                MCON(out T, out INFIN, out T, out T);
                CR = INFIN;
                CI = INFIN;
                return;

          _10:  if (Math.Abs(BR) >= Math.Abs(BI)) goto _20;
                R = BR / BI;
                D = BI + R * BR;
                CR = (AR * R + AI) / D;
                CI = (AI * R - AR) / D;
                return;

          _20:  R = BI / BR;
                D = BR + R * BI;
                CR = (AR + AI * R) / D;
                CI = (AI - AR * R) / D;
                return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A POLYNOMIAL.
        /// </summary>
        /// <param name="NN"></param>
        /// <param name="PT">PT IS THE MODULUS OF THE COEFFICIENTS.</param>
        /// <param name="Q"></param>
        /// <returns></returns>
        private double CAUCHY(int NN, double[] PT, double[] Q)
        {
            double FResult = double.NaN;
            try
            {
                double X, XM, F, DX, DF;
                PT[NN] = -PT[NN];

                // COMPUTE UPPER ESTIMATE OF BOUND.
                int N = NN - 1;             //?NN or NN-1
                X = Math.Exp( (Math.Log(-PT[NN]) - Math.Log(PT[C_FOR_START_IDX]) ) / (double)N ); //X = DEXP( (DLOG(-PT(NN)) - DLOG(PT(1)))/FLOAT(N) )
                if (PT[N] == 0.0) goto _20;

                // IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
                XM = -PT[NN] / PT[N];
                if (XM < X) X = XM;

                // CHOP THE INTERVAL (0,X) UNITL F LE 0.
          _20:  XM = X * 0.1;
                F = PT[C_FOR_START_IDX];
                for (int I = (C_FOR_START_IDX + 1); I <= NN; I++)    //DO 30 I = 2,NN
			    {
                    F = F * XM + PT[I];
			    }

                if (F <= 0.0) goto _40;
                if (double.IsNaN(F)) 
                    goto _40;

                X = XM;
                goto _20;
          _40:  DX = X;

                // DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES.
          _50:  if (Math.Abs(DX / X) <= 0.005) goto _70;

                Q[C_FOR_START_IDX] = PT[C_FOR_START_IDX];   //Q(1) = PT(1)

                for (int I = (C_FOR_START_IDX + 1); I <= NN; I++)    //DO 60 I = 2,NN
			    {
                    Q[I] = Q[I - 1] * X + PT[I];
			    }          
                  
                F = Q[NN];
                DF = Q[C_FOR_START_IDX];
                for (int I = (C_FOR_START_IDX + 1); I <= N; I++)     //DO 65 I = 2,N
			    {
                    DF = DF * X + Q[I];
			    }

                DX = F / DF;
                X = X - DX;

                //Avoid infinit loop.
                if (double.IsNaN(X))
                    goto _70;

                goto _50;

          _70:  FResult = X;
                return FResult;
            }
            catch (Exception ex)
            {
                throw ex;
            }

            //return FResult;
        }


        /// <summary>
        ///  COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
        ///  POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.
        /// </summary>
        /// <param name="L1">NUMBER OF NO-SHIFTS.</param>
        private void NOSHFT(int L1)
        {
            try
            {
                double XNI, T1, T2;
                int N, NM1, J;
                N = NN-1;
                NM1 = N-1;

                for (int I = C_FOR_START_IDX; I <= N; I++) //DO 10 I = 1,N
			    {
                    XNI = NN - I;
                    HR[I] = XNI * PR[I] / (double)N;
                    HI[I] = XNI * PI[I] / (double)N;
			    }

                for (int JJ = C_FOR_START_IDX; JJ <= L1; JJ++)     // DO 50 JJ = 1,L1
			    {
                    if (CMOD(HR[N], HI[N]) <= (ETA * 10.0 * CMOD(PR[N], PI[N]))) goto _30;

                    CDIVID(-PR[NN], -PI[NN], HR[N], HI[N], out TR, out TI);

                    for (int I = C_FOR_START_IDX; I <= NM1; I++)   //DO 20 I = 1,NM1
                    {
                        J = NN-I;
                        T1 = HR[J-1];
                        T2 = HI[J-1];
                        HR[J] = TR * T1-TI * T2 + PR[J];
                        HI[J] = TR * T2 + TI * T1 + PI[J];
                    }

                    HR[C_FOR_START_IDX] = PR[C_FOR_START_IDX];
                    HI[C_FOR_START_IDX] = PI[C_FOR_START_IDX];
                    goto _50;

                    // IF THE CONSTANT TERM IS ESSENTIALLY ZERO, SHIFT H COEFFICIENTS.
              _30:  for (int I = C_FOR_START_IDX; I <= NM1; I++) //DO 40 I = 1,NM1
			        {
                        J = NN-I;
                        HR[J] = HR[J - 1];
                        HI[J] = HI[J - 1];
			        }   
                    
                    HR[C_FOR_START_IDX] = 0.0;
                    HI[C_FOR_START_IDX] = 0.0;                    
			    }

          _50:  return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
        /// PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
        /// </summary>
        /// <param name="NN"></param>
        /// <param name="SR"></param>
        /// <param name="SI"></param>
        /// <param name="PR"></param>
        /// <param name="PI"></param>
        /// <param name="QR"></param>
        /// <param name="QI"></param>
        /// <param name="PVR"></param>
        /// <param name="PVI"></param>
        private void POLYEV(int NN, double SR, double SI, double[] PR, double[] PI, double[] QR, double[] QI, out double PVR, out double PVI)
        {
            try
            {
                double T;
                QR[C_FOR_START_IDX] = PR[C_FOR_START_IDX];
                QI[C_FOR_START_IDX] = PI[C_FOR_START_IDX];
                PVR = QR[C_FOR_START_IDX];
                PVI = QI[C_FOR_START_IDX];

                for (int I = (C_FOR_START_IDX + 1); I <= NN; I++)    //      DO 10 I = 2,NN
                {
                    T = PVR * SR - PVI * SI + PR[I];
                    PVI = PVR * SI + PVI * SR + PI[I];
                    PVR = T;
                    QR[I] = PVR;
                    QI[I] = PVI;
                }

                return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// COMPUTES  T = -P(S)/H(S).
        /// BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
        /// </summary>
        /// <param name="_BOOL">RETURNS TRUE IF H(S) IS ESSENTIALLY ZERO.</param>
        private void CALCT(out bool _BOOL)
        {
            try
            {
                double HVR, HVI;
                int N = NN - 1;

                // EVALUATE H(S).
                POLYEV(N, SR, SI, HR, HI, QHR, QHI, out HVR, out HVI);
                _BOOL = CMOD(HVR, HVI) <= (ARE * 10.0 * CMOD(HR[N], HI[N]));
                if (_BOOL) goto _10;

                CDIVID(-PVR, -PVI, HVR, HVI, out TR, out TI);
                return;

           _10: TR = 0.0;
                TI = 0.0;
                return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        ///  CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
        ///  _BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
        /// </summary>
        /// <param name="_BOOL">If true H[S] is essentially zero.</param>
        private void NEXTH(bool _BOOL)
        {
            try 
	        {
                double T1, T2;
                int N = NN - 1;
                //int NM1 = N - 1;
                if (_BOOL) goto _20;

                for (int J = (C_FOR_START_IDX + 1); J <= N; J++)     // DO 10 J = 2,N
                {
                    T1 = QHR[J - 1];
                    T2 = QHI[J - 1];
                    HR[J] = TR * T1 - TI * T2 + QPR[J];
                    HI[J] = TR * T2 + TI * T1 + QPI[J];
                }

                HR[C_FOR_START_IDX] = QPR[C_FOR_START_IDX];
                HI[C_FOR_START_IDX] = QPI[C_FOR_START_IDX];
                return;

                // IF H(S) IS ZERO REPLACE H WITH QH.
          _20:  for (int J = (C_FOR_START_IDX + 1); J <= N; J++)     //20 DO 30 J = 2,N
			    {
                     HR[J] = QHR[J - 1];
                     HI[J] = QHI[J - 1];
			    }
                   
                HR[C_FOR_START_IDX] = 0.0;
                HI[C_FOR_START_IDX] = 0.0;
                return;
	        }
	        catch (Exception ex)
	        {
		        throw ex;
	        }
        }


        /// <summary>
        /// CARRIES OUT THE THIRD STAGE ITERATION.
        /// </summary>
        /// <param name="L3">L3 - LIMIT OF STEPS IN STAGE 3.</param>
        /// <param name="ZR">ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.</param>
        /// <param name="ZI">ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.</param>
        /// <param name="CONV">RETURNS TRUE IF ITERATION CONVERGES.</param>
        private void VRSHFT(int L3, ref double ZR, ref double ZI, out bool CONV)
        {
            try
            {
                double MP, MS, OMP, RELSTP, R1, R2, TP;
                OMP = 0.0; RELSTP = 0.0;
                CONV = false;
                bool B = false;
                bool _BOOL;

                SR = ZR;
                SI = ZI;

                // MAIN LOOP FOR STAGE THREE
                for (int I = C_FOR_START_IDX; I <= L3; I++)    //      DO 60 I = 1,L3
                {
                    // EVALUATE P AT S AND TEST FOR CONVERGENCE.
                    POLYEV(NN, SR, SI, PR, PI, QPR, QPI, out PVR, out PVI);
                    MP = CMOD(PVR, PVI);
                    MS = CMOD(SR, SI);

                    if (MP > 20.0 * ERREV(NN, QPR, QPI, MS, MP, ARE, MRE)) goto _10;

                    // POLYNOMIAL VALUE IS SMALLER IN VALUE THAN A BOUND ON THE ERROR
                    // IN EVALUATING P, TERMINATE THE ITERATION.
                    CONV = true;
                    ZR = SR;
                    ZI = SI;
                    return;

              _10:  if (I == C_FOR_START_IDX) goto _40;   //TODO: if start at 0 put if (I==0)
                    if (B || (MP < OMP) || (RELSTP >= 0.05)) goto _30;
                         
                    // ITERATION HAS STALLED. PROBABLY A CLUSTER OF ZEROS. DO 5 FIXED
                    // SHIFT STEPS INTO THE CLUSTER TO FORCE ONE ZERO TO DOMINATE.
                    TP = RELSTP;
                    B = true;
                    
                    if (RELSTP < ETA) TP = ETA;

                    R1 = Math.Sqrt(TP);
                    R2 = SR * (1.0 + R1) - SI * R1;
                    SI = SR * R1 + SI * (1.0 + R1);
                    SR = R2;
                    POLYEV(NN, SR, SI, PR, PI, QPR, QPI, out PVR, out PVI);

                    for (int J = C_FOR_START_IDX; J <= 5; J++) //DO 20 J = 1,5
			        {
                        CALCT(out _BOOL);
                        NEXTH(_BOOL);
			        }

                    OMP = INFIN;
                    goto _50;
                    // EXIT IF POLYNOMIAL VALUE INCREASES SIGNIFICANTLY.
              _30:  if (MP * 0.1 > OMP) return;
              _40:  OMP = MP;
                    // CALCULATE NEXT ITERATE.
              _50:  CALCT(out _BOOL);
                    NEXTH(_BOOL);
                    CALCT(out _BOOL);
                    if (_BOOL) goto _60;
                    RELSTP = CMOD(TR, TI) / CMOD(SR, SI);
                    SR = SR + TR;
                    SI = SI + TI;
                }

          _60:  return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR
        /// CONVERGENCE.
        /// INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
        /// APPROXIMATE ZERO IF SUCCESSFUL.
        /// L2 - LIMIT OF FIXED SHIFT STEPS
        /// ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
        /// CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
        /// </summary>
        /// <param name="L2">Limit of fixed shift steps.</param>
        /// <param name="ZR">Return (if converge) zero real part.</param>
        /// <param name="ZI">Return (if converge) zero imaginary part.</param>
        /// <param name="CONV">Return indicating convergence or not of stage 3 iteration.</param>
        private void FXSHFT(int L2, ref double ZR, ref double ZI, out bool CONV)
        {
            try
            {
                double OTR, OTI, SVSR, SVSI;
                bool TEST, PASD, _BOOL;
                int N;

                N = NN - 1;

                // EVALUATE P AT S.
                POLYEV(NN, SR, SI, PR, PI, QPR, QPI, out PVR, out PVI);
                TEST = true;
                PASD = false;

                // CALCULATE FIRST T = -P(S)/H(S).
                CALCT(out _BOOL);

                //C MAIN LOOP FOR ONE SECOND STAGE STEP.
                for (int J = C_FOR_START_IDX; J <= L2; J++)    // DO 50 J = 1,L2
			    {
			        OTR = TR;
                    OTI = TI;

                    // COMPUTE NEXT H POLYNOMIAL AND NEW T.
                    NEXTH(_BOOL);
                    CALCT(out _BOOL);
                    ZR = SR + TR;
                    ZI = SI + TI;

                    // TEST FOR CONVERGENCE UNLESS STAGE 3 HAS FAILED ONCE OR THIS
                    // IS THE LAST H POLYNOMIAL .
                    if (_BOOL || !TEST || (J == L2)) goto _50;
                    if (CMOD(TR - OTR, TI - OTI) >= (0.5 * CMOD(ZR, ZI))) goto _40;
                    if (!PASD) goto _30;

                    // THE WEAK CONVERGENCE TEST HAS BEEN PASSED TWICE, START THE
                    // THIRD STAGE ITERATION, AFTER SAVING THE CURRENT H POLYNOMIAL
                    // AND SHIFT.
                    for (int I = C_FOR_START_IDX; I <= N; I++) //DO 10 I = 1,N;
			        {
			            SHR[I] = HR[I];
                        SHI[I] = HI[I];
			        }

                    SVSR = SR;
                    SVSI = SI;
                    VRSHFT(10, ref ZR, ref ZI, out CONV);
                    if (CONV) return;

                    // THE ITERATION FAILED TO CONVERGE. TURN OFF TESTING AND RESTORE
                    // H,S,PV AND T.
                    TEST = false;
                    for (int I = C_FOR_START_IDX; I <= N; I++) //DO 20 I = 1,N
			        {
                        HR[I] = SHR[I];
                        HI[I] = SHI[I];
			        }
                    
                    SR = SVSR;
                    SI = SVSI;
                    POLYEV(NN, SR, SI, PR, PI, QPR, QPI, out PVR, out PVI);
                    CALCT(out _BOOL);
                    goto _50;

              _30:  PASD = true;
                    goto _50;

              _40:  PASD = false;
			    }

              _50:  ;

                // ATTEMPT AN ITERATION WITH FINAL H POLYNOMIAL FROM SECOND STAGE.
                VRSHFT(10, ref ZR, ref ZI, out CONV);
                return;                          
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER
        /// RECURRENCE.
        /// </summary>
        /// <param name="NN"></param>
        /// <param name="QR">QR,QI - THE PARTIAL SUMS</param>
        /// <param name="QI">QR,QI - THE PARTIAL SUMS</param>
        /// <param name="MS">MS    -MODULUS OF THE POINT</param>
        /// <param name="MP">MP    -MODULUS OF POLYNOMIAL VALUE</param>
        /// <param name="ARE">ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION</param>
        /// <param name="MRE">ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION</param>
        /// <returns>RETURNS THE BOUNDS OF THE ERROR IN EVALUATING THE POLYNOMIAL.</returns>
        private double ERREV(int NN, double[] QR, double[] QI, double MS, double MP, double ARE, double MRE)
        {
            double FResult = double.NaN;
            try 
	        {	    
                double E;
		        E = CMOD(QR[C_FOR_START_IDX], QI[C_FOR_START_IDX]) * MRE / (ARE + MRE);     //      E = CMOD(QR(1),QI(1))*MRE/(ARE+MRE)

                for (int I = C_FOR_START_IDX; I <= NN; I++)    //      DO 10 I = 1,NN
			    {
                    E = E * MS + CMOD(QR[I], QI[I]);
			    }

                FResult = E * (ARE + MRE) - MP * MRE;               
	        }
	        catch (Exception ex)
	        {
		        throw ex;
	        }

            return FResult;
        }


        #endregion Private Methods



        /// <summary>
        /// Main routine to find polynomial zeros using Jenkin and Traub methods for polynomials with complex (or real) coeficents.
        /// </summary>
        /// <param name="OPR">  DOUBLE PRECISION VECTOR OF REAL PARTS OF THE COEFFICIENTS IN
        ///                     ORDER OF DECREASING POWERS.
        /// </param>
        /// <param name="OPI">  DOUBLE PRECISION VECTOR OF IMAGINARY PARTS OF THE COEFFICIENTS IN
        ///                     ORDER OF DECREASING POWERS.
        /// </param>
        /// <param name="DEGREE">INTEGER DEGREE OF POLYNOMIAL.</param>
        /// <param name="ZEROR">OUTPUT the real roots.</param>
        /// <param name="ZEROI">OUTPUT the imaginary roots.</param>
        /// <param name="FAIL"> OUTPUT LOGICAL PARAMETER,  TRUE  ONLY IF
        ///                     LEADING COEFFICIENT IS ZERO OR IF CPOLY
        ///                     HAS FOUND FEWER THAN DEGREE ZEROS.</param>
        /// <remarks>
        /// 
        /// IMPORTANT NOTE:
        ///     The coeficient arrays 'opr' (real part coeficients) and 'opi' (im part
        ///     coeficients) must have their first values at index 1 (not zero) because this
        ///     is a direct port from a FORTRAN77 algorithm and fortran uses one based array
        ///     indexing.
        /// 
        ///     THE PROGRAM HAS BEEN WRITTEN TO REDUCE THE CHANCE OF OVERFLOW
        ///     OCCURRING. IF IT DOES OCCUR, THERE IS STILL A POSSIBILITY THAT
        ///     THE ZEROFINDER WILL WORK PROVIDED THE OVERFLOWED QUANTITY IS
        ///     REPLACED BY A LARGE NUMBER.
        ///     
        /// </remarks>
        public void Roots(double[] OPR, double[] OPI, int DEGREE, out double[] ZEROR, out double[] ZEROI, out bool FAIL)
        {            
            try
            {
                int _size = DEGREE + 1; //Fortran arrays start with index 1 instead of 0.
                
                //Initialize out parameters zeros.
                FAIL = false;
                ZEROR = new double[_size];
                ZEROI = new double[_size];

                //Initialize globals
                PR = new double[_size + 1];
                PI = new double[_size + 1];
                HR = new double[_size];
                HI = new double[_size];
                QPR = new double[_size + 1];
                QPI = new double[_size + 1];
                QHR = new double[_size];
                QHI = new double[_size];
                SHR = new double[_size + 1];
                SHI = new double[_size + 1];

                // TO CHANGE THE SIZE OF POLYNOMIALS WHICH CAN BE SOLVED, REPLACE
                // THE DIMENSION OF THE ARRAYS IN THE COMMON AREA.

                double XX, YY, COSR, SINR, SMALNO, _BASE, XXX, ZR, ZI, BND;
                ZR = double.NaN;
                ZI = double.NaN;
                bool CONV = false;
                int CNT1, CNT2, IDNN2;

                //INITIALIZATION OF CONSTANTS
                this.MCON(out ETA, out INFIN, out SMALNO, out _BASE);

                ARE = ETA;
                MRE = 2.0 * Math.Sqrt(2.0) * ETA;// DSQRT(2.0D0)*ETA
                XX = 0.70710678;
                YY = -XX;
                COSR = -0.069756474;
                SINR = 0.99756405;
                FAIL = false;       // .FALSE.
                NN = DEGREE + 1;    // NN = DEGREE+1        TODO: CHECK THIS

                // ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO.
                if ((OPR[C_FOR_START_IDX] != 0.0) || (OPI[C_FOR_START_IDX] != 0.0))     //Arrays in FORTRAN start at 1 by default.
                    goto _10;

                FAIL = true;
                return;
              
                //REMOVE THE ZEROS AT THE ORIGIN IF ANY.
          _10:  if ((OPR[NN] != 0.0) || (OPI[NN] != 0.0)) 
                    goto _20;
                
                IDNN2 = DEGREE - NN + 2;
                ZEROR[IDNN2] = 0.0;     //Remmenber again that fortran arrays start at 1 by default.
                ZEROI[IDNN2] = 0.0;
                NN = NN - 1;
                goto _10;

                // MAKE A COPY OF THE COEFFICIENTS.
          _20:  for (int I = C_FOR_START_IDX; I <= NN; I++)   //20 DO 30 I = 1,NN
			    {
			        PR[I] = OPR[I];
                    PI[I] = OPI[I];
                    SHR[I] = CMOD(PR[I], PI[I]);
			    }

                // SCALE THE POLYNOMIAL.
                BND = SCALE(NN, SHR, ETA, INFIN, SMALNO,_BASE);
                if (BND == 1.0) goto _40;

                for (int I = C_FOR_START_IDX; I <= NN; I++)    //  DO 35 I = 1,NN
			    {
                    PR[I] = BND * PR[I];
                    PI[I] = BND * PI[I];
			    }

                // START THE ALGORITHM FOR ONE ZERO .
          _40:  if (NN > 2) goto _50;

                // CALCULATE THE FINAL ZERO AND RETURN.
                CDIVID(-PR[C_FOR_START_IDX + 1], -PI[C_FOR_START_IDX + 1], PR[C_FOR_START_IDX], PI[C_FOR_START_IDX], out ZEROR[DEGREE], out ZEROI[DEGREE]);
                return;
                
                // CALCULATE BND, A LOWER BOUND ON THE MODULUS OF THE ZEROS.
          _50:  for (int I = C_FOR_START_IDX; I <= NN; I++)   // 50 DO 60 I = 1,NN
			    {
                    SHR[I] = CMOD(PR[I], PI[I]);
			    }

                BND = CAUCHY(NN, SHR, SHI);

                // OUTER LOOP TO CONTROL 2 MAJOR PASSES WITH DIFFERENT SEQUENCES
                // OF SHIFTS.
                for (CNT1 = 1; CNT1 <= 2; CNT1++)    //DO 100 CNT1 = 1,2
			    {
                    // FIRST STAGE CALCULATION, NO SHIFT.
                    NOSHFT(5);

                    // INNER LOOP TO SELECT A SHIFT.
                    for (CNT2 = 1; CNT2 <= 9; CNT2++)    // DO 90 CNT2 = 1,9
                    {
                        // SHIFT IS CHOSEN WITH MODULUS BND AND AMPLITUDE ROTATED BY
                        // 94 DEGREES FROM THE PREVIOUS SHIFT
                        XXX = COSR * XX - SINR * YY;
                        YY = SINR * XX + COSR * YY;
                        XX = XXX;
                        SR = BND * XX;
                        SI = BND * YY;

                        // SECOND STAGE CALCULATION, FIXED SHIFT.
                        FXSHFT(10 * CNT2, ref ZR, ref ZI, out CONV);
                        if (!CONV) goto _80;

                        // THE SECOND STAGE JUMPS DIRECTLY TO THE THIRD STAGE ITERATION.
                        // IF SUCCESSFUL THE ZERO IS STORED AND THE POLYNOMIAL DEFLATED.
                        IDNN2 = DEGREE - NN + 2;
                        ZEROR[IDNN2] = ZR;
                        ZEROI[IDNN2] = ZI;
                        NN = NN - 1;

                        for (int I = C_FOR_START_IDX; I <= NN; I++)    //DO 70 I = 1,NN
			            {
                            PR[I] = QPR[I];
                            PI[I] = QPI[I];
			            }
                        
                        goto _40;
                  _80:  ;
                        // IF THE ITERATION IS UNSUCCESSFUL ANOTHER SHIFT IS CHOSEN.
                     
                    }
                    
                    // IF 9 SHIFTS FAIL, THE OUTER LOOP IS REPEATED WITH ANOTHER
                    // SEQUENCE OF SHIFTS.                    
			    }
                
                // THE ZEROFINDER HAS FAILED ON TWO MAJOR PASSES.
                // RETURN EMPTY HANDED.
                FAIL = true;
                return;
            }
            catch (Exception)
            {
                throw;
            }
        }



    }


}

