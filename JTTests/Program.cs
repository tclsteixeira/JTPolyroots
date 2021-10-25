using System;

namespace JTTests
{
    class MainClass
    {

        public static void Main(string[] args)
        {
            Console.WriteLine("--------------------------------");

            Console.WriteLine("TESTING JENKINS-TRAUB ALGORITHMS FOR FINDING POLYNOMIAL ROOTS");
            Console.WriteLine("RPOLY - FOR POLYNOMIALS WITH REAL COEFICIENTS ONLY");
            Console.WriteLine("CPOLY - FOR POLYNOMIALS WITH REAL AND COMPLEX COEFICIENTS");

            Console.WriteLine("--------------------------------");
            Console.WriteLine();

            TestRPoly.Test();

            Console.WriteLine();
            Console.WriteLine();

            TestCPoly.Test();

        }

    }
}
