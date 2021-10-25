using System;

namespace JTPolyroots
{

    /// <summary>
    /// Jenkins-Traub real (double precision) coeficients polynomial root finder.
    /// --- Applies only to one variable polynomials with real (double precision) coeficients ---
    /// 
    /// </summary>
    /// <remarks>
    /// This is a direct port from FORTRAN90.
    /// It works for polynomials only with Real coeficents.
    /// If you are using only polynomials with Real coeficients use instead RPoly because it is 4 times faster than CPoly.
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
    public sealed class JT_RPoly
    {

        //Start index for fortran arrays.
        private const int C_FOR_START_IDX = 1;



        #region globals

        //double[] p, qp, k, qk, svk, temp, pt;
        double[] p, qp, k, qk, svk;
        double eta, szr, szi, lzr, lzi, are, mre;
        double sr, si, u, v;

        double a, b, c, d, e, f, g, h;
        double a1, a7, a3;

        int n, nn;

        #endregion globals



        public JT_RPoly()
            : base()
        { }



        #region Private methods


        /// <summary>
        /// Calculate the zeros of the quadratic a*z**2+b1*z+c.
        /// The quadratic formula, modified to avoid overflow, is used to find the
        /// larger zero if the zeros are real and both zeros are complex.
        /// The smaller real zero is found directly from the product of the zeros c/a.
        /// </summary>
        /// <param name="a">a parameter.</param>
        /// <param name="b1">b parameter.</param>
        /// <param name="c">c parameter.</param>
        /// <param name="sr"></param>
        /// <param name="si"></param>
        /// <param name="lr"></param>
        /// <param name="li"></param>
        private void quad(double a, double b1, double c, out double sr, out double si, out double lr, out double li)
        {
            try
            {
                double b, d, e;
                if (a != 0.0) goto _20;

                sr = 0.0;
                if (b1 != 0.0) sr = -c / b1;

                lr = 0.0;
          _10:  si = 0.0;
                li = 0.0;
                return;

          _20:  if (c == 0.0) 
                {
                  sr = 0.0;
                  lr = -b1 / a;
                  goto _10;
                }

                //! Compute discriminant avoiding overflow
                b = b1 / 2.0;
                if (Math.Abs(b) >= Math.Abs(c)) 
                {
                    e = 1.0 - (a/b) * (c/b);
                    d = Math.Sqrt(Math.Abs(e)) * Math.Abs(b);
                }
                else
                {
                    e = a;
                    if (c < 0.0) e = -a;

                    e = b * (b / Math.Abs(c)) - e;
                    d = Math.Sqrt(Math.Abs(e)) * Math.Sqrt(Math.Abs(c));
                }

                if (e >= 0.0) 
                {
                    // Real zeros
                    if (b >= 0.0) d = -d;

                    lr = (-b + d) / a;
                    sr = 0.0;
                    if (lr != 0.0) sr = (c / lr) / a;

                    goto _10;
                }
                
                //! complex conjugate zeros
                sr = -b / a;
                lr = sr;
                si = Math.Abs(d / a);
                li = -si;
                return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Computes up to  l2  fixed shift k-polynomials, testing for convergence in
        /// the linear or quadratic case.  Initiates one of the variable shift
        /// iterations and returns with the number of zeros found.
        /// l2 - limit of fixed shift steps
        /// nz - number of zeros found
        /// </summary>
        /// <param name="l2">Limit of fixed shift steps.</param>
        /// <param name="nz">Returns the number of zeros found.</param>
        private void fxshfr(int l2, out int nz)
        {
            try 
	        {
                double svu, svv, ui, vi, s;
                double betas, betav, oss, ovv, ss, vv, ts, tv, ots, otv, tvv, tss;
                otv = 0.0; ots = 0.0;
                int _TYPE, j, iflag;
                bool vpass, spass, vtry, stry;

                nz = 0;
                betav = 0.25;
                betas = 0.25;
                oss = sr;
                ovv = v;

                // Evaluate polynomial by synthetic division
                quadsd(nn, u, v, p, out qp, out a, out b);
                calcsc(out _TYPE);
                for (j = C_FOR_START_IDX; j <= l2; j++)    //DO  j = 1, l2
                {
                    // calculate next k polynomial and estimate v
                    nextk(_TYPE);
                    calcsc(out _TYPE);
                    newest(_TYPE, out ui, out vi);
                    vv = vi;

                    // Estimate s
                    ss = 0.0;
                    if (k[n] != 0.0) ss = -p[nn] / k[n];

                    tv = 1.0;
                    ts = 1.0;
                    if ((j != 1) && (_TYPE != 3))
                    {
                        // Compute relative measures of convergence of s and v sequences
                        if (vv != 0.0) tv = Math.Abs((vv - ovv) / vv);

                        if (ss != 0.0) ts = Math.Abs((ss - oss) / ss);

                        // If decreasing, multiply two most recent convergence measures
                        tvv = 1.0;
                        if (tv < otv) tvv = tv * otv;

                        tss = 1.0;
                        if (ts < ots) tss = ts * ots;

                        // Compare with convergence criteria
                        vpass = tvv < betav;
                        spass = tss < betas;
                        if (spass || vpass)
                        {
                            // At least one sequence has passed the convergence test.
                            // Store variables before iterating
                            svu = u;
                            svv = v;
                            CopyArray(k, svk, C_FOR_START_IDX, n);  // svk(1:n) = k(1:n)                            
                            s = ss;

                            // Choose iteration according to the fastest converging sequence
                            vtry = false;
                            stry = false;
                            if (spass && ((!vpass) || (tss < tvv))) goto _40;

                        _20: quadit(ui, vi, out nz);
                            if (nz > 0) return;

                            // Quadratic iteration has failed. flag that it has
                            // been tried and decrease the convergence criterion.
                            vtry = true;
                            betav = betav * 0.25;

                            // Try linear iteration if it has not been tried and
                            // the s sequence is converging
                            if (stry || (!spass)) goto _50;

                            CopyArray(svk, k, C_FOR_START_IDX, n);    //  k(1:n) = svk(1:n)
                      _40:  realit(ref s, out nz, out iflag);
                            if (nz > 0) return;

                            // Linear iteration has failed.  Flag that it has been
                            // tried and decrease the convergence criterion
                            stry = true;
                            betas = betas * 0.25;
                            if (iflag != 0)
                            {
                                // If linear iteration signals an almost double real
                                // zero attempt quadratic interation
                                ui = -(s + s);
                                vi = s * s;
                                goto _20;
                            }

                            // Restore variables
                        _50: u = svu;
                            v = svv;
                            CopyArray(svk, k, C_FOR_START_IDX, n);  // k(1:n) = svk(1:n)

                            // Try quadratic iteration if it has not been tried
                            // and the v sequence is converging
                            if (vpass && (!vtry)) goto _20;

                            // Recompute qp and scalar values to continue the second stage
                            quadsd(nn, u, v, p, out qp, out a, out b);
                            calcsc(out _TYPE);
                        }
                    }

                    ovv = vv;
                    oss = ss;
                    otv = tv;
                    ots = ts;
                }

                return;
	        }
	        catch (Exception ex)
	        {
		        throw ex;
	        }
        }


        
        /// <summary>
        /// Divides p by the quadratic  1,u,v  placing the
        /// quotient in q and the remainder in a,b.
        /// </summary>
        /// <param name="?"></param>
        /// <param name="?"></param>
        /// <param name="?"></param>
        /// <param name="?"></param>
        /// <param name="?"></param>
        /// <param name="?"></param>
        /// <param name="?"></param>
        private void quadsd(int nn, double u, double v, double[] p, out double[] q, out double a, out double b)
        {
            try 
	        {
                double c;
                int i;
                q = new double[nn + C_FOR_START_IDX];

                b = p[C_FOR_START_IDX];
                q[C_FOR_START_IDX] = b;
                a = p[C_FOR_START_IDX + 1] - u * b;
                q[C_FOR_START_IDX + 1] = a;
                for (i = (C_FOR_START_IDX + 2); i <= nn; i++)   //DO  i = 3, nn
                {
                    c = p[i] - u * a - v * b;
                    q[i] = c;
                    b = a;
                    a = c;
                }

                return;
	        }
	        catch (Exception ex)
	        {
		        throw ex;
	        }
        }


        /// <summary>
        /// This routine calculates scalar quantities used to
        /// compute the next k polynomial and new estimates of
        /// the quadratic coefficients.
        /// type - integer variable set here indicating how the
        /// calculations are normalized to avoid overflow
        /// </summary>
        /// <param name="_TYPE">
        /// integer variable set here indicating how the
        /// calculations are normalized to avoid overflow.
        /// </param>
        private void calcsc(out int _TYPE)
        {
            try
            {
                /// Synthetic division of k by the quadratic 1,u,v
                quadsd(n, u, v, k, out qk, out c, out d);
                if (Math.Abs(c) <= Math.Abs(k[n]) * 100.0 * eta)
                {
                    if (Math.Abs(d) <= Math.Abs(k[n - 1]) * 100.0 * eta)
                    {
                        _TYPE = 3;
                        // type=3 indicates the quadratic is almost a factor of k
                        return;
                    }
                }

                if (Math.Abs(d) >= Math.Abs(c))
                {
                    _TYPE = 2;

                    // type=2 indicates that all formulas are divided by d
                    e = a / d;
                    f = c / d;
                    g = u * b;
                    h = v * b;
                    a3 = (a + g) * e + h * (b / d);
                    a1 = b * f - a;
                    a7 = (f + u) * a + h;
                    return;
                }

                _TYPE = 1;
                // type=1 indicates that all formulas are divided by c
                e = a / c;
                f = d / c;
                g = u * e;
                h = v * b;
                a3 = a * e + (h/c+g) * b;
                a1 = b - a * (d/c);
                a7 = a + g * d + h * f;
                return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Computes the next k polynomials using scalars computed in calcsc.
        /// </summary>
        /// <param name="_TYPE"></param>
        private void nextk(int _TYPE)
        {
            try
            {
                double temp;
                int i;
                
                if (_TYPE != 3) 
                {
                    temp = a;
                    if (_TYPE == 1) temp = b;

                    if (Math.Abs(a1) <= Math.Abs(temp)*eta*10.0) 
                    {
                        // If a1 is nearly zero then use a special form of the recurrence
                        k[C_FOR_START_IDX] = 0.0;
                        k[C_FOR_START_IDX + 1] = -a7 * qp[C_FOR_START_IDX];
                        for (i = (C_FOR_START_IDX + 2); i <= n; i++)    // DO  i = 3, n
                        {
                            k[i] = a3 * qk[i-2] - a7 * qp[i-1];
                        }

                        return;
                    }

                    //! Use scaled form of the recurrence
                    a7 = a7 / a1;
                    a3 = a3 / a1;
                    k[C_FOR_START_IDX] = qp[C_FOR_START_IDX];
                    k[C_FOR_START_IDX + 1] = qp[C_FOR_START_IDX + 1] - a7 * qp[C_FOR_START_IDX];
                    for (i = (C_FOR_START_IDX + 2); i <= n; i++)  //DO  i = 3, n
                    {
                        k[i] = a3 * qk[i - 2] - a7 * qp[i - 1] + qp[i];
                    }
                    return;
                }

                // Use unscaled form of the recurrence if type is 3
                k[C_FOR_START_IDX] = 0.0;
                k[C_FOR_START_IDX + 1] = 0.0;
                for (i = (C_FOR_START_IDX + 2); i <= n; i++)  //DO  i = 3, n
                {
                    k[i] = qk[i - 2];
                }
                return;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Compute new estimates of the quadratic coefficients
        /// using the scalars computed in calcsc.
        /// </summary>
        /// <param name="_TYPE"></param>
        /// <param name="uu"></param>
        /// <param name="vv"></param>
        private void newest(int _TYPE, out double uu, out double vv)
        {
            try 
	        {
                double a4, a5, b1, b2, c1, c2, c3, c4, temp;

                // Use formulas appropriate to setting of type.
                if (_TYPE != 3)
                {
                    if (_TYPE != 2)
                    {
                        a4 = a + u * b + h * f;
                        a5 = c + (u + v * f) * d;
                    }
                    else
                    {
                        a4 = (a + g) * f + h;
                        a5 = (f + u) * c + v * d;
                    }
                
                    // Evaluate new quadratic coefficients.
                    b1 = -k[n] / p[nn];
                    b2 = -(k[n-1] + b1 * p[n]) / p[nn];
                    c1 = v * b2 * a1;
                    c2 = b1 * a7;
                    c3 = b1 * b1 * a3;
                    c4 = c1 - c2 - c3;
                    temp = a5 + b1 * a4 - c4;
                    if (temp != 0.0) 
                    {
                        uu = u - (u*(c3+c2)+v*(b1*a1+b2*a7)) / temp;
                        vv = v * (1.0 + c4/temp);
                        return;
                    }
                }

                // If type=3 the quadratic is zeroed
                uu = 0.0;
                vv = 0.0;
                return;
	        }
	        catch (Exception ex)
	        {
		        throw ex;
	        }
        }


        /// <summary>
        /// Variable-shift k-polynomial iteration for a quadratic factor, converges
        /// only if the zeros are equimodular or nearly so.
        /// uu,vv - coefficients of starting quadratic
        /// nz - number of zero found
        /// </summary>
        /// <param name="uu">Coefficients of starting quadratic.</param>
        /// <param name="vv">Coefficients of starting quadratic</param>
        /// <param name="nz">Number of zero found.</param>
        private void quadit(double uu, double vv, out int nz)
        {
            try
            {
                double ui, vi;
                double mp, omp, ee, relstp, t, zm;
                int _TYPE, i, j;
                bool tried;
                relstp = 0.0;
                omp = 0.0;

                nz = 0;
                tried = false;
                u = uu;
                v = vv;
                j = 0;

                // Main loop
          _10:  quad(1.0, u, v, out szr, out szi, out lzr, out lzi);

                // Return if roots of the quadratic are real and not
                // close to multiple or nearly equal and  of opposite sign.
                if (Math.Abs(Math.Abs(szr) - Math.Abs(lzr)) > 0.01 * Math.Abs(lzr)) return;

                // Evaluate polynomial by quadratic synthetic division
                quadsd(nn, u, v, p, out qp, out a, out b);
                mp = Math.Abs(a - szr * b) + Math.Abs(szi * b);

                // Compute a rigorous  bound on the rounding error in evaluting p
                zm = Math.Sqrt(Math.Abs( (double)v) );
                ee = 2.0 * Math.Abs( (double)(qp[C_FOR_START_IDX]) );
                t = -szr * b;

                for (i = (C_FOR_START_IDX + 1); i <= n; i++)    //DO  i = 2, n
                {
                    ee = ee * zm + Math.Abs( (double)(qp[i]) );
                }

                ee = ee * zm + Math.Abs( (double)(a) + t );
                ee = (5.0 * mre + 4.0 * are) * ee - (5.0 * mre + 2.0 * are) * (Math.Abs( (double)(a) + t) +
                     Math.Abs( (double)(b) ) * zm) + 2.0 * are * Math.Abs(t);

                // Iteration has converged sufficiently if the
                // polynomial value is less than 20 times this bound
                if (mp <= 20.0 * ee)
                {
                    nz = 2;
                    return;
                }

                j = j + 1;

                // Stop iteration after 20 steps
                if (j > 20) return;
                if (j >= 2)
                {
                    if (!((relstp > 0.01) || (mp < omp) || tried))
                    {
                        // A cluster appears to be stalling the convergence.
                        // five fixed shift steps are taken with a u,v close to the cluster
                        if (relstp < eta) relstp = eta;

                        relstp = Math.Sqrt(relstp);
                        u = u - u * relstp;
                        v = v + v * relstp;
                        quadsd(nn, u, v, p, out qp, out a, out b);
                        for (i = 1; i <= 5; i++)//    DO  i = 1, 5
                        {
                            calcsc(out _TYPE);
                            nextk(_TYPE);
                        }

                        tried = true;
                        j = 0;
                    }
                }

                omp = mp;

                //! Calculate next k polynomial and new u and v
                calcsc(out _TYPE);
                nextk(_TYPE);
                calcsc(out _TYPE);
                newest(_TYPE, out ui, out vi);

                // If vi is zero the iteration is not converging
                if (vi == 0.0) return;

                relstp = Math.Abs((vi-v)/vi);
                u = ui;
                v = vi;
                goto _10;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Variable-shift h polynomial iteration for a real zero.
        /// sss   - starting iterate
        /// nz    - number of zero found
        /// iflag - flag to indicate a pair of zeros near real axis.
        /// </summary>
        /// <param name="sss">Starting iterate.</param>
        /// <param name="nz">Number of zero found.</param>
        /// <param name="iflag">Flag to indicate a pair of zeros near real axis.</param>
        private void realit(ref double sss, out int nz, out int iflag)
        {
            try 
	        {
                double pv, kv, t, s;
                t = 0.0;
                double  ms, mp, omp, ee;
                omp = 0.0;
                int i, j;

                nz = 0;
                s = sss;
                iflag = 0;
                j = 0;

                //! Main loop
          _10:  pv = p[C_FOR_START_IDX];

                //! Evaluate p at s
                qp[C_FOR_START_IDX] = pv;
                for (i = (C_FOR_START_IDX + 1); i <= nn; i++)   //DO  i = 2, nn
                {
                    pv = pv * s + p[i];
                    qp[i] = pv;
                }

                mp = Math.Abs(pv);

                // Compute a rigorous bound on the error in evaluating p
                ms = Math.Abs(s);
                ee = (mre/(are+mre)) * Math.Abs( (double)(qp[C_FOR_START_IDX]));
                for (i = (C_FOR_START_IDX + 1); i <= nn; i++)  //DO  i = 2, nn
                {
                    ee = ee * ms + Math.Abs((double)(qp[i]));
                }

                // Iteration has converged sufficiently if the
                // polynomial value is less than 20 times this bound
                if (mp <= 20.0 * ((are+mre)*ee - mre*mp)) 
                {
                    nz = 1;
                    szr = s;
                    szi = 0.0;
                    return;
                }

                j = j + 1;

                // Stop iteration after 10 steps
                if (j > 10) return;

                if (j >= 2)
                {
                    if (Math.Abs(t) <= 0.001 * Math.Abs(s - t) && mp > omp)
                    {
                        // A cluster of zeros near the real axis has been encountered,
                        // return with iflag set to initiate a quadratic iteration
                        iflag = 1;
                        sss = s;
                        return;
                    }
                }

                // Return if the polynomial value has increased significantly
                omp = mp;

                // Compute t, the next polynomial, and the new iterate
                kv = k[C_FOR_START_IDX];
                qk[C_FOR_START_IDX] = kv;
                for (i = (C_FOR_START_IDX + 1); i <= n; i++)    //DO  i = 2, n
                {
                    kv = kv * s + k[i];
                    qk[i] = kv;
                }

                if (Math.Abs(kv) > (Math.Abs(k[n]) * 10.0 * eta)) 
                {
                    // Use the scaled form of the recurrence if the value of k at s is nonzero
                    t = -pv / kv;
                    k[C_FOR_START_IDX] = qp[C_FOR_START_IDX];
                    for (i = (C_FOR_START_IDX + 1); i <= n; i++)  //DO  i = 2, n
                    {
                        k[i] = t * qk[i-1] + qp[i];
                    }
                }
                else
                {
                    //! Use unscaled form
                    k[C_FOR_START_IDX] = 0.0;
                    for (i = (C_FOR_START_IDX + 1); i <= n ; i++)  //  DO  i = 2, n
                    {
                        k[i] = qk[i-1];
                    }
                }

                kv = k[C_FOR_START_IDX];
                for (i=(C_FOR_START_IDX + 1); i <= n; i++)  //DO  i = 2, n
                {
                    kv = kv * s + k[i];
                }

                t = 0.0;
                if (Math.Abs(kv) > (Math.Abs(k[n]) * 10.0 * eta)) t = -pv / kv;

                s = s + t;
                goto _10;
	        }
	        catch (Exception ex)
	        {
		        throw ex;
	        }
        }
        

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
                INFINY = _BASE * (1.0E0 - Math.Pow(_BASE, -T)) * Math.Pow(_BASE, M - 1);
                SMALNO = Math.Pow(_BASE, N + 3) / Math.Pow(_BASE, 3);
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }


        /// <summary>
        /// Copies the source array elements to a destination array starting at a given index until a given last index.
        /// </summary>
        /// <param name="source">The source array.</param>
        /// <param name="dest">The destination array.</param>
        /// <param name="start">The start index of source where copy begins.</param>
        /// <param name="last">The last index of source where copy stops (included).</param>
        private void CopyArray(double[] source, double[] dest, int start, int last)
        {
            try
            {
                Array.Copy(source, start, dest, start, last - start + 1);
            }
            catch (Exception ex)
            {
                throw ex;
            }            
        }


        #endregion Private methods




        /// <summary>
        ///! Finds the zeros of a real polynomial
        ///! op  - double precision vector of coefficients in order of
        ///!       decreasing powers.
        ///! degree   - integer degree of polynomial.
        ///! zeror, zeroi - output double precision vectors of real and imaginary parts
        ///!                of the zeros.
        ///! fail  - output logical parameter, true only if leading coefficient is zero
        ///!         or if rpoly has found fewer than degree zeros.
        ///!         In the latter case degree is reset to the number of zeros found.
        ///! To change the size of polynomials which can be solved, reset the dimensions
        ///! of the arrays in the common area and in the following declarations.
        ///! The subroutine uses single precision calculations for scaling, bounds and
        ///! error calculations.  All calculations for the iterations are done in
        ///! double precision.        /// 
        /// </summary>
        /// <param name="op">double precision vector of coefficients in order of decreasing powers.</param>
        /// <param name="degree">integer degree of polynomial.</param>
        /// <param name="zeror">output double precision vectors of real parts of the zeros.</param>
        /// <param name="zeroi">output double precision vectors of imaginary parts of the zeros.</param>
        /// <param name="fail">output logical parameter, true only if leading coefficient is zero
        ///                    or if rpoly has found fewer than degree zeros.
        ///                    In the latter case degree is reset to the number of zeros found.
        /// </param>
        /// <returns>The number of roots found.</returns>
        public int Roots(double[] op, int degree, out double[] zeror, out double[] zeroi, out bool fail)
        {
            int FResult = degree;
            try
            {

                //! COMMON /global/ p, qp, k, qk, svk, sr, si, u, v, a, b, c, d, a1,  &
                //!    a2, a3, a6, a7, e, f, g, h, szr, szi, lzr, lzi, eta, are, mre, n, nn

                //double[] p, qp, k, qk, svk;
                //double sr, si, u, v, a, b, c, d, a1, a2, a3, a6, 
                //       a7, e, f, g, h, szr, szi, lzr, lzi; 


                double lo, MAX, MIN, xx, yy, cosr, sinr, xxx, x, sc, bnd,
               xm, ff, df, dx, infin, smalno, _base;
                //eta, szr, szi, lzr, lzi;
                double[] temp, pt;
                double l;   //don´t know if shoud be integer??
                int nm1, cnt, nz, j, i, jj;
                double factor, aa, bb, cc, t;
                bool zerok;
                int _size = degree + 1; //Fortran arrays start with index 1 instead of 0.

                //Initialize out parameters zeros.
                //FAIL = false;
                zeror = new double[_size];
                zeroi = new double[_size];

                //INITIALIZATION OF CONSTANTS
                this.MCON(out eta, out infin, out smalno, out _base);

                //! are and mre refer to the unit error in + and * respectively.
                //! They are assumed to be the same as eta.
                are = eta;
                mre = eta;
                lo = smalno / eta;

                //! Initialization of constants for shift rotation
                xx = Math.Sqrt(0.5);
                yy = -xx;
                cosr = -.069756474;
                sinr = .99756405;
                fail = false;
                n = degree;
                nn = n + 1;

                //! Algorithm fails if the leading coefficient is zero.
                if (op[1] == 0.0)
                {
                    fail = true;
                    degree = 0;
                    FResult = degree;
                    return FResult;
                }

            //! Remove the zeros at the origin if any
            _10: if (op[nn] == 0.0)
                {
                    j = degree - n + 1;
                    zeror[j] = 0.0;
                    zeroi[j] = 0.0;
                    nn = nn - 1;
                    n = n - 1;
                    goto _10;
                }

                //! Allocate various arrays
                p = new double[nn + 1]; qp = new double[nn + 1]; k = new double[nn + 1]; qk = new double[nn + 1];
                svk = new double[nn + 1]; temp = new double[nn + 1]; pt = new double[nn + 1];

                //! Make a copy of the coefficients
                this.CopyArray(op, p, C_FOR_START_IDX, nn);     // p[1:nn] = op[1:nn]

            //! Start the algorithm for one zero
            _30: if (n <= 2)
                {
                    if (n < 1) return FResult;

                    //! calculate the final zero or pair of zeros
                    if (n != 2)
                    {
                        zeror[degree] = -p[2] / p[1];
                        zeroi[degree] = 0.0;
                        FResult = degree;
                        return FResult;
                    }

                    quad(p[1], p[2], p[3], out zeror[degree - 1], out zeroi[degree - 1],
                         out zeror[degree], out zeroi[degree]);
                    FResult = degree;
                    return FResult;
                }

                //! Find largest and smallest moduli of coefficients.
                MAX = 0.0;
                MIN = infin;

                for (i = C_FOR_START_IDX; i <= nn; i++) //DO  i = 1, nn
                {
                    x = Math.Abs((double)(p[i]));
                    if (x > MAX) MAX = x;
                    if (x != 0.0 && x < MIN) MIN = x;
                }

                //! Scale if there are large or very small coefficients computes a scale
                //! factor to multiply the coefficients of the polynomial.
                //! The scaling is done to avoid overflow and to avoid undetected underflow
                //! interfering with the convergence criterion.
                //! The factor is a power of the base
                sc = lo / MIN;

                if (sc <= 1.0)
                {
                    if (MAX < 10.0) goto _60;
                    if (sc == 0.0) sc = smalno;
                }
                else
                {
                    if ((infin / sc) < MAX) goto _60;
                }

                l = Math.Log(sc) / Math.Log(_base) + 0.5;
                factor = Math.Pow((_base * 1.0), l);

                if (factor != 1.0)
                {
                    for (i = C_FOR_START_IDX; i <= nn; i++)  //    p[1:nn) = factor * p(1:nn)
                    {
                        p[i] = factor * p[i];
                    }
                }

            //! compute lower bound on moduli of zeros.

            _60: for (i = C_FOR_START_IDX; i <= nn; i++) //60 pt(1:nn) = ABS(p(1:nn))
                {
                    pt[i] = Math.Abs(p[i]);
                }

                pt[nn] = -pt[nn];

                //! compute upper estimate of bound
                x = Math.Exp((Math.Log(-pt[nn]) - Math.Log(pt[C_FOR_START_IDX])) / n);
                if (pt[n] != 0.0)
                {
                    //! if newton step at the origin is better, use it.
                    xm = -pt[nn] / pt[n];
                    if (xm < x) x = xm;
                }

            //! chop the interval (0,x) until ff .le. 0
            _80: xm = x * 0.1;
                ff = pt[1];
                for (i = (C_FOR_START_IDX + 1); i <= nn; i++)    // DO  i = 2, nn
                {
                    ff = ff * xm + pt[i];
                }

                if (ff > 0.0)
                {
                    x = xm;
                    goto _80;
                }

                dx = x;

            //! do newton iteration until x converges to two decimal places
            _100: if (Math.Abs(dx / x) > 0.005)
                {
                    ff = pt[C_FOR_START_IDX];
                    df = ff;

                    for (i = (C_FOR_START_IDX + 1); i <= n; i++)    //DO  i = 2, n
                    {
                        ff = ff * x + pt[i];
                        df = df * x + ff;
                    }

                    ff = ff * x + pt[nn];
                    dx = ff / df;
                    x = x - dx;
                    goto _100;
                }

                bnd = x;

                //! compute the derivative as the intial k polynomial
                //! and do 5 steps with no shift
                nm1 = n - 1;
                for (i = (C_FOR_START_IDX + 1); i <= n; i++)    //DO  i = 2, n
                {
                    k[i] = (nn - i) * p[i] / n;
                }

                k[C_FOR_START_IDX] = p[C_FOR_START_IDX];
                aa = p[nn];
                bb = p[n];
                zerok = (k[n] == 0.0);

                for (jj = 1; jj <= 5; jj++) //DO  jj = 1, 5
                {
                    cc = k[n];
                    if (!zerok)
                    {
                        //! use scaled form of recurrence if value of k at 0 is nonzero
                        t = -aa / cc;
                        for (i = C_FOR_START_IDX; i <= nm1; i++)    //    DO  i = 1, nm1
                        {
                            j = nn - i;
                            k[j] = t * k[j - 1] + p[j];
                        }

                        k[C_FOR_START_IDX] = p[C_FOR_START_IDX];
                        zerok = Math.Abs(k[n]) <= Math.Abs(bb) * eta * 10.0;
                    }
                    else
                    {
                        //! use unscaled form of recurrence
                        for (i = C_FOR_START_IDX; i <= nm1; i++)  //    DO  i = 1, nm1
                        {
                            j = nn - i;
                            k[j] = k[j - 1];
                        }

                        k[C_FOR_START_IDX] = 0.0;
                        zerok = (k[n] == 0.0);
                    }
                }

                //! save k for restarts with new shifts
                for (i = C_FOR_START_IDX; i <= n; i++)  // temp(1:n) = k(1:n)
                {
                    temp[i] = k[i];
                }

                //! loop to select the quadratic  corresponding to each
                //! new shift
                for (cnt = C_FOR_START_IDX; cnt <= 20; cnt++) //DO  cnt = 1, 20
                {
                    //! Quadratic corresponds to a double shift to a non-real point and its complex
                    //! conjugate.  The point has modulus bnd and amplitude rotated by 94 degrees
                    //! from the previous shift
                    xxx = cosr * xx - sinr * yy;
                    yy = sinr * xx + cosr * yy;
                    xx = xxx;
                    sr = bnd * xx;
                    si = bnd * yy;
                    u = -2.0 * sr;
                    v = bnd;

                    //! second stage calculation, fixed quadratic
                    fxshfr(20 * cnt, out nz);
                    if (nz != 0)
                    {
                        //! The second stage jumps directly to one of the third stage iterations and
                        //! returns here if successful.
                        //! Deflate the polynomial, store the zero or zeros and return to the main
                        //! algorithm.
                        j = degree - n + 1;
                        zeror[j] = szr;
                        zeroi[j] = szi;
                        nn = nn - nz;
                        n = nn - 1;
                        CopyArray(qp, p, C_FOR_START_IDX, nn);  //p(1:nn) = qp(1:nn);

                        if (nz == 1) goto _30;
                        zeror[j + 1] = lzr;
                        zeroi[j + 1] = lzi;
                        goto _30;
                    }

                    //! If the iteration is unsuccessful another quadratic
                    //! is chosen after restoring k
                    CopyArray(temp, k, C_FOR_START_IDX, nn);   //k(1:nn) = temp(1:nn)
                }

                //! Return with failure if no convergence with 20 shifts
                fail = true;
                degree = degree - n;
                FResult = degree;
                return FResult;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }



    }


}

