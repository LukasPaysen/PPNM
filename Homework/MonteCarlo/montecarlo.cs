// MonteCarlo.cs
// Plain Monte Carlo integration + Halton quasi-Monte Carlo for the singular integral
// Produces data files for plotting error scaling.

using System;
using System.IO;
using System.Collections.Generic;

class MonteCarlo
{
    // single Random instance for pseudo-random sampling
    static Random rnd = new Random(12345);

    // Plain pseudo-random Monte Carlo as in the exercise
    static (double, double) plainmc(Func<double[], double> f, double[] a, double[] b, int N)
    {
        int dim = a.Length;
        double V = 1.0;
        for (int i = 0; i < dim; i++) V *= (b[i] - a[i]);

        double sum = 0.0, sum2 = 0.0;
        double[] x = new double[dim];

        for (int i = 0; i < N; i++)
        {
            for (int k = 0; k < dim; k++) x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
            double fx = f(x);
            sum += fx;
            sum2 += fx * fx;
        }

        double mean = sum / N;
        double varianceTerm = sum2 / N - mean * mean;
        if (varianceTerm < 0 && varianceTerm > -1e-15) varianceTerm = 0; // guard tiny negative
        double sigma = Math.Sqrt(Math.Max(0.0, varianceTerm));
        double integral = mean * V;
        double error = sigma * V / Math.Sqrt(N);
        return (integral, error);
    }

    // Halton low-discrepancy sequence (radical inverse in different prime bases)
    static readonly int[] primes = new int[] {2,3,5,7,11,13,17,19,23,29,31};

    static double RadicalInverse(int index, int @base)
    {
        // index is >= 1 (use 1-based to avoid returning 0 for index=0)
        double result = 0.0;
        double f = 1.0 / @base;
        int i = index;
        while (i > 0)
        {
            int digit = i % @base;
            result += digit * f;
            i /= @base;
            f /= @base;
        }
        return result;
    }

    static double[] HaltonPoint(int index, int dim)
    {
        double[] u = new double[dim];
        for (int k = 0; k < dim; k++)
        {
            int b = primes[k];
            u[k] = RadicalInverse(index, b);
        }
        return u;
    }

    // Quasi-Monte Carlo integration using Halton sequence starting at given index
    // We return the plain estimate (integral) and also the sample variance if wanted.
    static (double, double) quasiPlainmc(Func<double[], double> f, double[] a, double[] b, int N, int startIndex)
    {
        int dim = a.Length;
        double V = 1.0;
        for (int i = 0; i < dim; i++) V *= (b[i] - a[i]);

        double sum = 0.0, sum2 = 0.0;
        double[] x = new double[dim];

        // Halton indices: use startIndex + i (1-based index advised)
        for (int i = 0; i < N; i++)
        {
            double[] u = HaltonPoint(startIndex + i + 1, dim); // +1 so index 0 -> 1
            for (int k = 0; k < dim; k++) x[k] = a[k] + u[k] * (b[k] - a[k]);
            double fx = f(x);
            sum += fx;
            sum2 += fx * fx;
        }

        double mean = sum / N;
        double varianceTerm = sum2 / N - mean * mean;
        if (varianceTerm < 0 && varianceTerm > -1e-15) varianceTerm = 0;
        double sigma = Math.Sqrt(Math.Max(0.0, varianceTerm));
        double integral = mean * V;
        double pseudoError = sigma * V / Math.Sqrt(N); // not a statistically correct error for QMC, but informative
        return (integral, pseudoError);
    }

    // Useful test integrands
    static double CircleIndicator(double[] x)
    {
        double xx = x[0], yy = x[1];
        return (xx * xx + yy * yy <= 1.0) ? 1.0 : 0.0;
    }

    static double PolyFunc(double[] x)
    {
        return x[0] * x[0] + x[1] * x[1];
    }

    // The difficult singular integrand
    // Note: the assignment uses measure dx/pi etc. We'll integrate over [0,pi]^3 and divide by pi^3
    static double SingularFunc(double[] x)
    {
        double cx = Math.Cos(x[0]);
        double cy = Math.Cos(x[1]);
        double cz = Math.Cos(x[2]);
        double denom = 1.0 - cx * cy * cz;
        // guard against (very) small denom -- return a large number (integrand is integrable but peaked)
        if (denom <= 0.0) return 1e300; // avoid crash, but in practice denom>0 on (0,pi]
        return 1.0 / denom;
    }

    static void Main(string[] args)
    {
        // --- Part A tests (2D) ---
        int[] Ns2D = new int[] { 100, 500, 1000, 5000, 10000, 50000, 100000 };

        string circleFile = "circle.dat";
        string polyFile = "poly.dat";

        // Unit circle area over [-1,1]x[-1,1]
        using (var sc = new StreamWriter(circleFile))
        {
            sc.WriteLine("# N	estimate	estimated_error	actual_error	1/sqrt(N)");
            double[] a = new double[] { -1.0, -1.0 };
            double[] b = new double[] { 1.0, 1.0 };
            double exactCircle = Math.PI;
            foreach (int N in Ns2D)
            {
                var res = plainmc(CircleIndicator, a, b, N);
                double est = res.Item1;
                double estErr = res.Item2;
                double actErr = Math.Abs(est - exactCircle);
                sc.WriteLine($"{N}	{est:F12}	{estErr:E6}	{actErr:E6}	{1.0 / Math.Sqrt(N):E6}");
                Console.WriteLine($"Circle N={N}: est={est:F12}, estErr={estErr:E6}, actErr={actErr:E6}");
            }
        }

        // Polynomial integral over [0,1]x[0,1]
        using (var sp = new StreamWriter(polyFile))
        {
            sp.WriteLine("# N	estimate	estimated_error	actual_error	1/sqrt(N)");
            double[] a2 = new double[] { 0.0, 0.0 };
            double[] b2 = new double[] { 1.0, 1.0 };
            double exactPoly = 2.0 / 3.0; // integral of x^2+y^2 over unit square
            foreach (int N in Ns2D)
            {
                var res = plainmc(PolyFunc, a2, b2, N);
                double est = res.Item1;
                double estErr = res.Item2;
                double actErr = Math.Abs(est - exactPoly);
                sp.WriteLine($"{N}	{est:F12}	{estErr:E6}	{actErr:E6}	{1.0 / Math.Sqrt(N):E6}");
                Console.WriteLine($"Poly N={N}: est={est:F12}, estErr={estErr:E6}, actErr={actErr:E6}");
            }
        }

        // --- Part B: Singular integral experiments (3D) ---
        int[] Ns3D = new int[] { 1000, 5000, 10000, 50000, 100000 };

        string singularPlainFile = "singular_plain.dat";
        string singularQuasiFile = "singular_quasi.dat";
        string singularCompareFile = "singular_compare.dat";

        // exact value (given) for reference
        double exactSingular = 1.3932039296856768591842462603255; // high-precision constant

        // Domain for the integral
        double pi = Math.PI;
        double[] a3 = new double[] { 0.0, 0.0, 0.0 };
        double[] b3 = new double[] { pi, pi, pi };

        // Plain Monte Carlo on the singular integral
        using (var sp = new StreamWriter(singularPlainFile))
        {
            sp.WriteLine("# N	estimate	estimated_error	actual_error	1/sqrt(N)");
            foreach (int N in Ns3D)
            {
                var res = plainmc(SingularFunc, a3, b3, N);
                double integral = res.Item1; // this approximates integral over [0,pi]^3 of F(x) dx
                // desired integral includes factor (1/pi)^3 -> divide by pi^3, but as noted mean = integral/V,
                // and integral/(pi^3) = mean, so simply compute:
                double estimate = integral / (pi * pi * pi);
                double estErr = res.Item2 / (pi * pi * pi);
                double actErr = Math.Abs(estimate - exactSingular);
                sp.WriteLine($"{N}	{estimate:F12}	{estErr:E6}	{actErr:E6}	{1.0/Math.Sqrt(N):E6}");
                Console.WriteLine($"Singular plain N={N}: est={estimate:F12}, estErr={estErr:E6}, actErr={actErr:E6}");
            }
        }

        // Quasi-Monte Carlo (Halton) on the singular integral
        // We'll generate two Halton-based estimates with different start indices and use their difference as an error proxy
        using (var sq = new StreamWriter(singularQuasiFile))
        using (var scmp = new StreamWriter(singularCompareFile))
        {
            sq.WriteLine("# N	estimate1	estimate2	diff_estimate	actual_err1	actual_err2	1/sqrt(N)");
            scmp.WriteLine("# N	plain_act_err	quasi_act_err	plain_est_err	quasi_est_err	1/sqrt(N)");

            int offset1 = 1; // start index for first Halton sequence
            int offset2 = 100000; // far offset for second sequence to get a different quasi-sample

            foreach (int N in Ns3D)
            {
                var q1 = quasiPlainmc(SingularFunc, a3, b3, N, offset1);
                var q2 = quasiPlainmc(SingularFunc, a3, b3, N, offset2);

                double integral1 = q1.Item1 / (pi * pi * pi);
                double integral2 = q2.Item1 / (pi * pi * pi);

                double estDiff = Math.Abs(integral1 - integral2); // quasi estimated error proxy
                double actErr1 = Math.Abs(integral1 - exactSingular);
                double actErr2 = Math.Abs(integral2 - exactSingular);

                sq.WriteLine($"{N}	{integral1:F12}	{integral2:F12}	{estDiff:E6}	{actErr1:E6}	{actErr2:E6}	{1.0/Math.Sqrt(N):E6}");
                scmp.WriteLine($"{N}	{GetPlainActErrFromFile(singularPlainFile,N):E6}	{Math.Min(actErr1,actErr2):E6}	{GetPlainEstErrFromFile(singularPlainFile,N):E6}	{estDiff:E6}	{1.0/Math.Sqrt(N):E6}");

                Console.WriteLine($"Singular quasi N={N}: est1={integral1:F12}, est2={integral2:F12}, diff={estDiff:E6}, actErr(min)={Math.Min(actErr1,actErr2):E6}");
            }
        }

        Console.WriteLine("
Data written to: circle.dat, poly.dat, singular_plain.dat, singular_quasi.dat, singular_compare.dat");
        Console.WriteLine("Use the plotting script to visualize actual vs estimated error and compare with 1/sqrt(N).");
    }

    // Small helpers to read back plainfile entries for comparison (simple and robust)
    static double GetPlainActErrFromFile(string filename, int N)
    {
        if (!File.Exists(filename)) return double.NaN;
        foreach (var line in File.ReadLines(filename))
        {
            if (line.StartsWith("#")) continue;
            var parts = line.Split(new char[] {'	',' ', ','}, StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length < 4) continue;
            if (int.TryParse(parts[0], out int n0) && n0 == N)
            {
                if (double.TryParse(parts[3], out double val)) return val;
            }
        }
        return double.NaN;
    }

    static double GetPlainEstErrFromFile(string filename, int N)
    {
        if (!File.Exists(filename)) return double.NaN;
        foreach (var line in File.ReadLines(filename))
        {
            if (line.StartsWith("#")) continue;
            var parts = line.Split(new char[] {'	',' ', ','}, StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length < 3) continue;
            if (int.TryParse(parts[0], out int n0) && n0 == N)
            {
                if (double.TryParse(parts[2], out double val)) return val;
            }
        }
        return double.NaN;
    }
}

/*
Makefile (put in a file named 'Makefile')

run:
	mcs -out:MonteCarlo.exe MonteCarlo.cs
	mono MonteCarlo.exe

clean:
	rm -f MonteCarlo.exe circle.dat poly.dat singular_plain.dat singular_quasi.dat singular_compare.dat

*/
