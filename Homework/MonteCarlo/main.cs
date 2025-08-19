// MonteCarlo.cs
// Plain Monte Carlo integration + Halton quasi-Monte Carlo for the singular integral
// Produces data files for plotting error scaling.

using System;
using System.IO;
using System.Collections.Generic;

class MonteCarlo
{
    // ===== Stratified sampling (N-point, recursive) =====

// Base: use our existing plainmc
static (double val, double err) stratified_mc(
    Func<double[], double> f, double[] a, double[] b, int N, int nmin)
{
    int dim = a.Length;
    if (N <= nmin) return plainmc(f, a, b, N);

    // --- pilot sampling (nmin points) to decide where to split ---
    double V = 1.0; for (int k = 0; k < dim; k++) V *= (b[k] - a[k]);
    int m = nmin;

    double sum = 0.0;
    int[] nL = new int[dim], nR = new int[dim];
    double[] sumL = new double[dim], sumR = new double[dim];
    double[] x = new double[dim];

    for (int i = 0; i < m; i++)
    {
        for (int k = 0; k < dim; k++) x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
        double fx = f(x);
        sum += fx;
        for (int k = 0; k < dim; k++)
        {
            double midk = 0.5 * (a[k] + b[k]);
            if (x[k] >= midk) { nR[k]++; sumR[k] += fx; }
            else              { nL[k]++; sumL[k] += fx; }
        }
    }

    double mean = sum / m;
    // sub-variance proxy per dimension: |mean_right - mean_left|
    int kdiv = 0; double best = -1.0;
    for (int k = 0; k < dim; k++)
    {
        double mL = (nL[k] > 0) ? (sumL[k] / nL[k]) : mean;
        double mR = (nR[k] > 0) ? (sumR[k] / nR[k]) : mean;
        double varProxy = Math.Abs(mR - mL);
        if (varProxy > best) { best = varProxy; kdiv = k; }
    }

    // Make child boxes split at midpoint along kdiv
    double[] aL = (double[])a.Clone(), bL = (double[])b.Clone();
    double[] aR = (double[])a.Clone(), bR = (double[])b.Clone();
    double mid = 0.5 * (a[kdiv] + b[kdiv]);
    bL[kdiv] = mid; aR[kdiv] = mid;

    // Estimate sub-variances for allocation (use same proxy but recomputed for chosen dim)
    double mLdiv = (nL[kdiv] > 0) ? (sumL[kdiv] / nL[kdiv]) : mean;
    double mRdiv = (nR[kdiv] > 0) ? (sumR[kdiv] / nR[kdiv]) : mean;
    double vL = Math.Abs(mLdiv - mean) + 1e-16;
    double vR = Math.Abs(mRdiv - mean) + 1e-16;

    // Allocate remaining points (we already spent m). Ensure both sides get at least nmin/2 when possible.
    int Nrem = Math.Max(0, N - m);
    int Nleft = (int)Math.Round(Nrem * (vL / (vL + vR)));
    int Nright = Nrem - Nleft;

    // Guard: if either side got too few and base wouldn’t trigger, bump minimally
    int minChild = Math.Max(2, nmin / 2);
    if (Nleft > 0 && Nleft < minChild && Nright > minChild) { Nleft = minChild; Nright = Nrem - Nleft; }
    if (Nright > 0 && Nright < minChild && Nleft > minChild) { Nright = minChild; Nleft = Nrem - Nright; }

    // Recurse; add the pilot m points implicitly by trusting child estimates (standard approach)
    var left  = (Nleft > 0)  ? stratified_mc(f, aL, bL, Nleft, nmin)  : (0.0, 0.0);
    var right = (Nright > 0) ? stratified_mc(f, aR, bR, Nright, nmin) : (0.0, 0.0);

    // Combine: integral is additive; errors combine in quadrature (assume independence)
    double val = left.Item1 + right.Item1;
    double err = Math.Sqrt(left.Item2 * left.Item2 + right.Item2 * right.Item2);

    // If allocation gave all to one side (degenerate case), fall back to plain on the whole box
    if (Nrem > 0 && Nleft == 0 && Nright == 0)
        return plainmc(f, a, b, N);

    return (val, err);
}

    // single Random instance for pseudo-random sampling
    static Random rnd = new Random(100);
    // --- add these helpers to MonteCarlo class ---

// HaltonPointShifted: returns Halton point (index) shifted by 'shift' (wrap mod 1)
static double[] HaltonPointShifted(int index, int dim, double[] shift)
{
    double[] u = HaltonPoint(index, dim);
    for (int k = 0; k < dim; k++)
    {
        u[k] = u[k] + shift[k];
        if (u[k] >= 1.0) u[k] -= 1.0;
    }
    return u;
}

// quasiWithShift: single Halton run using provided shift vector; returns integral estimate (not divided by pi^3)
static (double,double) quasiWithShift(Func<double[],double> f, double[] a, double[] b, int N, int startIndex, double[] shift)
{
    int dim = a.Length;
    double V = 1.0;
    for (int i = 0; i < dim; i++) V *= (b[i] - a[i]);

    double sum = 0.0, sum2 = 0.0;
    double[] x = new double[dim];

    for (int i = 0; i < N; i++)
    {
        double[] u = HaltonPointShifted(startIndex + i + 1, dim, shift);
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
    double pseudoErr = sigma * V / Math.Sqrt(N); // just informative
    return (integral, pseudoErr);
}

// RepeatQuasiShift: do R independent random shifts, return (meanEstimate, stddevEstimate)
// 'seed' controls reproducibility (use fixed for repeatable runs, or Time-dependent for varied runs)
static (double mean, double stddev) RepeatQuasiShift(Func<double[],double> f, double[] a, double[] b, int N, int startIndex, int R, int seed = 5555)
{
    int dim = a.Length;
    double[] estimates = new double[R];
    var rndLocal = new Random(seed);

    for (int r = 0; r < R; r++)
    {
        double[] shift = new double[dim];
        for (int k = 0; k < dim; k++) shift[k] = rndLocal.NextDouble();
        var res = quasiWithShift(f, a, b, N, startIndex, shift);
        estimates[r] = res.Item1;
    }

    double m = 0; foreach (var v in estimates) m += v; m /= R;
    double s = 0; foreach (var v in estimates) s += (v - m) * (v - m); s = Math.Sqrt(s / Math.Max(1, R - 1));
    return (m, s);
}


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
        Console.WriteLine("Task A: \nWe calculate two-dimensional integrals with our Monte-Carlo routine. We calculate the area of the unit circle and the poly integral ∬_{[0,1]×[0,1]} (x^2 + y^2) dx dy.\n");

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

        Console.WriteLine("\nThe estimated error and the actual error as functions of the number of sampling points along with 1/sqrt(N) have been plotted in the files circle_errors.png and poly_errors.png for the unit area and the poly integral respectively.\n");
        Console.WriteLine("Now we try to calculate the singular integral with our pseudo-random Monte-Carlo routine:\n");

        // --- Part B: Singular integral experiments (3D) ---
        int[] Ns3D = new int[] {100, 1000, 5000, 10000, 50000, 100000 };

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

        Console.WriteLine("\nTask B:");
        Console.WriteLine("We now Implement a multidimensional Monte-Carlo integrator that uses low-discrepancy (quasi-random) sequences. \n");

        // Example usage (inside Main) replacing the earlier quasi block:
    int R = 12;
    int startIndex = 1;
    int seedForShifts = 5555;

    using (var sq = new StreamWriter("singular_quasi_rand.dat"))
    using (var scmp = new StreamWriter("singular_compare_rand.dat"))
    {
        sq.WriteLine("# N\\tmean_estimate\\tstddev_estimate\\tactual_err_mean\\t1/sqrt(N)");
        scmp.WriteLine("# N\\tplain_act_err\\tquasi_act_err_mean\\tplain_est_err\\tquasi_stddev\\t1/sqrt(N)");

        foreach (int N in Ns3D)
        {
            // Repeat quasi with random shifts
            var qres = RepeatQuasiShift(SingularFunc, a3, b3, N, startIndex, R, seedForShifts);
            double meanIntegral = qres.mean / (pi * pi * pi);
            double stddev = qres.stddev / (pi * pi * pi);

            double plainActErr = GetPlainActErrFromFile("singular_plain.dat", N); // unchanged
            double plainEstErr = GetPlainEstErrFromFile("singular_plain.dat", N); // unchanged
            double quasiActErr = Math.Abs(meanIntegral - exactSingular);

            sq.WriteLine($"{N}\t{meanIntegral:F12}\t{stddev:E6}\t{quasiActErr:E6}\t{1.0/Math.Sqrt(N):E6}");
            scmp.WriteLine($"{N}\t{plainActErr:E6}\t{quasiActErr:E6}\t{plainEstErr:E6}\t{stddev:E6}\t{1.0/Math.Sqrt(N):E6}");

            Console.WriteLine($"Quasi-shifted N={N}: mean={meanIntegral:F12}, stddev={stddev:E6}, actErr={quasiActErr:E6}");
        }
    }


        Console.WriteLine("\nWe Compare the scaling of the error with our pseudo-random Monte-Carlo integrator in the file singular_compare_rand.png \nAlthough the quasi-random method does not clearly outperform the pseudo-random in terms of accuracy, we see that Plain Monte Carlo is very noisy for this singular integral, with large outliers when samples hit spike regions. In contrast, random-shifted Halton QMC produces smoother, more stable estimates with smaller standard deviations across shifts.");
        
        Console.WriteLine("\nTask C: \nWe compare stratified sampling to our original Monte-Carlo routine:\n");
        // --- Part C: Stratified vs Plain on 2D unit circle ---
int[] NsStrata = new int[] { 200, 500, 1000, 5000, 10000, 50000 };
int nmin = 64; // small pilot size (tweak 32..128 if you like)

using (var sf = new StreamWriter("strata_compare_circle.dat"))
{
    sf.WriteLine("# N\tplain_est\tplain_acterr\tstrata_est\tstrata_err\tstrata_acterr\t1/sqrt(N)");
    double[] aC = new double[] { -1.0, -1.0 };
    double[] bC = new double[] {  1.0,  1.0 };
    double exact = Math.PI;

    foreach (int N in NsStrata)
    {
        var pr = plainmc(CircleIndicator, aC, bC, N);
        double pest = pr.Item1, pact = Math.Abs(pest - exact);

        var sr = stratified_mc(CircleIndicator, aC, bC, N, nmin);
        double sest = sr.Item1, sact = Math.Abs(sest - exact);

        sf.WriteLine($"{N}\t{pest:F12}\t{pact:E6}\t{sest:F12}\t{sr.Item2:E6}\t{sact:E6}\t{1.0/Math.Sqrt(N):E6}");
        Console.WriteLine($"Strata circle N={N}: plain actErr={pact:E3}, strata actErr={sact:E3}, strata estErr={sr.Item2:E3}");
    }
}
Console.WriteLine("\nThe plot where we compare stratified sampling vs plain MC can be seen in the file strata_circle.png. Here the trend is even more obvious - the plain MC is not necessarily more inaccurate but is more noizy with large outliers.");

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
