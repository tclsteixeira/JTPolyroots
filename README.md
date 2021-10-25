## JTPolyroots  
  
  
This is the home page of **JTPolyroots**, a C# port from fortran libraries [CPOLY](https://calgo.acm.org/419.gz) and [RPOLY](https://calgo.acm.org/493.gz).  
**CPOLY** and **RPOLY** are self-contained numeric libraries that provides an efficient and accurate implementation of [Jenkins and Traub algorithm to compute polynomial roots](https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm). 
They gave two variants, one for general polynomials with complex coefficients, commonly known as the "CPOLY" algorithm, and a more complicated variant for the special case of polynomials with real coefficients, commonly known as the "RPOLY" algorithm. The latter is "practically a standard in black-box polynomial root-finders".
It's one of the most valuable method for polynomial root finding used for ex. in [Maple](https://www.maplesoft.com/) and [Mathematica](https://www.wolfram.com/mathematica/) software.

CPOLY is implemented in the Roots(...) procedure of the JT_CPoly class, while RPoly is implemented in the Roots(...) procedure of the JT_RPoly class.

### Usage example

    using system;
    using JTPolyroots;
    
    namespace Test
    {
        class  MainClass  
        {
            public  static  void Main(string[]  args)  
            {  
                // Note: on polynomial coeficients arrays, coeficients start at index one 
                //       and must be in descending order
                Console.WriteLine("--------------------------------");  
                JT_CPoly  cpoly  =  new  JT_CPoly();  
                try  
                {  
                    int degree = 3;
                    int size = degree + 2; // Note: first index start at 1
                    double[]  _P  =  new  double[size];  // coeficients real part
                    double[]  _PI  =  new  double[size];  // coeficients imaginary part
                    double[]  _ZR  =  new  double[size];  // roots real part
                    double[]  _ZI  =  new  double[size];  // roots imaginary part
                    bool  _FAIL  =  false;  
                    string  msgFail = "CPoly has failed on this example";  
                    System.Console.WriteLine("Solving third degree polynomial 'x^3-55x^2+1320x-18150'");  
                    
                    // Set polynomial coeficients (real part)
                    _P[1]  =  1;  
                    _P[2]  =  -55;  
                    _P[3]  =  1320;  
                    _P[4]  =  -18150;  
            
                    // test with imaginary coeficient values equal to zero
                    for  (int  I  =  1;  I  <  size;  I++) 
                    {  
                        _PI[I]  =  0.0;  
                    }  
  
                    // compute roots
                    cpoly.Roots(_P,  _PI,  degree,  out  _ZR,  out  _ZI,  out  _FAIL);

                    if  (_FAIL)  
                    {  
                        System.Console.WriteLine(msgFail);  
                    }
                    else
                    {   
                        System.Console.WriteLine("Roots:");
                        // print results
                        for (int I = 1; I < size; I++)
                        {
                            System.Console.WriteLine("{0} + ({1})i", _ZR[I], _ZI[I])
                        }
                    }           
            
                    System.Console.WriteLine();
                    System.Console.WriteLine("Press <Enter> to exit..."); 
                    System.Console.ReadLine();
                }
            }
        }
    }
  
    
**Author JTPolyroots:**  
<i>Copyright (c) 2021 Tiago C. Teixeira</i>  
  
**Algorithm copyright:**

<i>Copyright (c) [Michael A. Jenkins](https://research.cs.queensu.ca/home/maj/)</i>
<i>Copyright (c) [Joseph F. Traub](http://dli.library.cmu.edu/traub/)</i>

## Accuracy  
  
The methods have been extensively tested by many people. As predicted they enjoy faster than quadratic convergence for all distributions of zeros.

However, there are polynomials which can cause loss of precision as illustrated by the following example. The polynomial has all its zeros lying on two half-circles of different radii. [Wilkinson](https://en.wikipedia.org/wiki/James_H._Wilkinson "James H. Wilkinson") recommends that it is desirable for stable deflation that smaller zeros be computed first. The second-stage shifts are chosen so that the zeros on the smaller half circle are found first. After deflation the polynomial with the zeros on the half circle is known to be ill-conditioned if the degree is large; see Wilkinson,[[9]](https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm#cite_note-9) p. 64. The original polynomial was of degree 60 and suffered severe deflation instability.
  
## Copyright and Citation from original Fortran library  

<i>Copyright © 2021 [Association for Computing Machinery, Inc.](https://www.acm.org/)  </i>

Permission to include in application software or to make digital or hard copies of part  
or all of this work is subject to the following licensing agreement:
[ACM Software License Agreement](https://www.acm.org/publications/policies/software-copyright-notice)

When using CPoly and RPoly in scientific work, please cite as follows:  

* ALGORITHM 419 COLLECTED ALGORITHMS FROM ACM.  
* ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 02,  
P. 097.  
 

Copyright (C) [Michael A. Jenkins](https://research.cs.queensu.ca/home/maj/), Queen's University, Kingston, Canada, 1972, 1975; [Joseph F. Traub](http://dli.library.cmu.edu/traub/), Carnegie Mellon University, 1972, 1975.  
  
License: [ACM Software License Agreement](https://www.acm.org/publications/policies/software-copyright-notice)
  

## Further references  
  
Added changes from Remark on Algorithm 419 by David H. Withers  
CACM (March 1974) Vol 17 No 3 p. 157
  
  
> Please report bugs to the package maintainer.  
  
> Written with [StackEdit](https://stackedit.io/).
