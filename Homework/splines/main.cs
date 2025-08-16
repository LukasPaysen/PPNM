using System;
using System.IO;
using System.Globalization;

class Program
{
    // Provided binary search (exactly as in the assignment)
    public static int binsearch(double[] x, double z)
    {
        /* locates the interval for z by bisection */
        if (z < x[0] || z > x[x.Length - 1]) throw new Exception("binsearch: bad z");
        int i = 0, j = x.Length - 1;
        while (j - i > 1)
        {
            int mid = (i + j) / 2;
            if (z > x[mid]) i = mid; else j = mid;
        }
        return i;
    }

    // Linear interpolation using binsearch
    public static double linterp(double[] x, double[] y, double z)
    {
        int i = binsearch(x, z);
        double dx = x[i + 1] - x[i]; if (!(dx > 0)) throw new Exception("uups...");
        double dy = y[i + 1] - y[i];
        return y[i] + dy / dx * (z - x[i]);
    }

    // Analytical integral of the linear spline from x[0] to z.
    // We sum full trapezoid areas for intervals before the one containing z
    // and add the partial trapezoid for the final interval.
    public static double linterpInteg(double[] x, double[] y, double z)
    {
        if (z < x[0] || z > x[x.Length - 1]) throw new Exception("linterpInteg: bad z");
        int i = binsearch(x, z);
        double integral = 0.0;
        // full intervals
        for (int k = 0; k < i; ++k)
        {
            double dx = x[k + 1] - x[k];
            // area of trapezoid = (y[k] + y[k+1]) / 2 * dx
            integral += 0.5 * (y[k] + y[k + 1]) * dx;
        }
        // partial interval [x[i], z]
        double dxp = z - x[i];
        if (dxp > 0)
        {
            double y_at_z = linterp(x, y, z);
            integral += 0.5 * (y[i] + y_at_z) * dxp;
        }
        return integral;
    }

    public class Qspline
    {
        double[] x, y, b, c;
        double[] prefIntegral;
        int n;

        public Qspline(double[] xs, double[] ys)
        {
            if (xs.Length != ys.Length) throw new Exception("xs/ys length mismatch");
            n = xs.Length;
            if (n < 2) throw new Exception("need at least two points");
            x = (double[])xs.Clone();
            y = (double[])ys.Clone();
            int m = n - 1;
            b = new double[m];
            c = new double[m];

            double[] h = new double[m];
            double[] d = new double[m];
            double[] slope = new double[m];
            for (int i = 0; i < m; ++i)
            {
                h[i] = x[i + 1] - x[i];
                if (!(h[i] > 0)) throw new Exception("xs must be strictly increasing");
                d[i] = y[i + 1] - y[i];
                slope[i] = d[i] / h[i];
            }

            // Boundary: c0 = 0 => b0 = slope0
            b[0] = slope[0];
            // recurrence b[i+1] = 2*slope[i] - b[i]
            for (int i = 0; i < m - 1; ++i) b[i + 1] = 2.0 * slope[i] - b[i];

            // compute c[i]
            for (int i = 0; i < m; ++i) c[i] = (d[i] - b[i] * h[i]) / (h[i] * h[i]);

            // prefix integrals: prefIntegral[k] = integral from x0 to x[k]
            prefIntegral = new double[n];
            prefIntegral[0] = 0.0;
            for (int i = 0; i < m; ++i)
            {
                double H = h[i];
                double full = y[i] * H + 0.5 * b[i] * H * H + (1.0 / 3.0) * c[i] * H * H * H;
                prefIntegral[i + 1] = prefIntegral[i] + full;
            }
        }

        // evaluate spline at z
        public double Evaluate(double z)
        {
            int i = Program.binsearch(x, z); // uses your existing binsearch
            double t = z - x[i];
            return y[i] + b[i] * t + c[i] * t * t;
        }

        // derivative
        public double Derivative(double z)
        {
            int i = Program.binsearch(x, z);
            double t = z - x[i];
            return b[i] + 2.0 * c[i] * t;
        }

        // definite integral from x[0] to z
        public double Integral(double z)
        {
            int i = Program.binsearch(x, z);
            double t = z - x[i];
            double contrib = y[i] * t + 0.5 * b[i] * t * t + (1.0 / 3.0) * c[i] * t * t * t;
            return prefIntegral[i] + contrib;
        }

        // debug accessors
        public double[] GetB() { return (double[])b.Clone(); }
        public double[] GetC() { return (double[])c.Clone(); }
    }

    public class CubicSpline
    {
        double[] x, y, M;          // node arrays and second derivatives M[i]=S''(x_i)
        double[] A, B, C, D;      // per-interval coefficients for S_i(t)=A + B t + C t^2 + D t^3
        double[] prefIntegral;    // prefix integrals: prefIntegral[k] = integral from x0 to x[k]
        int n;

        public CubicSpline(double[] xs, double[] ys)
        {
            if (xs.Length != ys.Length) throw new Exception("xs/ys length mismatch");
            n = xs.Length;
            if (n < 2) throw new Exception("need at least two points");
            x = (double[])xs.Clone();
            y = (double[])ys.Clone();

            int m = n - 1;
            double[] h = new double[m];
            for (int i = 0; i < m; ++i)
            {
                h[i] = x[i + 1] - x[i];
                if (!(h[i] > 0)) throw new Exception("xs must be strictly increasing");
            }

            if (n == 2)
            {
                // degenerate: single interval -> linear / trivial cubic with M=0
                M = new double[n];
                M[0] = M[1] = 0.0;
            }
            else
            {
                // build tridiagonal system for interior M[1..n-2]
                int msys = n - 2;
                double[] a = new double[msys]; // sub-diagonal
                double[] diag = new double[msys];
                double[] c = new double[msys]; // super-diagonal
                double[] rhs = new double[msys];

                for (int i = 0; i < msys; ++i)
                {
                    int ii = i + 1; // corresponds to node index ii
                    a[i] = h[ii - 1];
                    diag[i] = 2.0 * (h[ii - 1] + h[ii]);
                    c[i] = h[ii];
                    double s1 = (y[ii + 1] - y[ii]) / h[ii];
                    double s0 = (y[ii] - y[ii - 1]) / h[ii - 1];
                    rhs[i] = 6.0 * (s1 - s0);
                }

                // Thomas algorithm (compact)
                double[] cp = new double[msys];
                double[] dp = new double[msys];
                cp[0] = c[0] / diag[0];
                dp[0] = rhs[0] / diag[0];
                for (int i = 1; i < msys; ++i)
                {
                    double denom = diag[i] - a[i] * cp[i - 1];
                    cp[i] = (i == msys - 1) ? 0.0 : c[i] / denom; // cp[last] not used but set to 0
                    dp[i] = (rhs[i] - a[i] * dp[i - 1]) / denom;
                }
                double[] Mint = new double[msys];
                Mint[msys - 1] = dp[msys - 1];
                for (int i = msys - 2; i >= 0; --i) Mint[i] = dp[i] - cp[i] * Mint[i + 1];

                // assemble full M with natural BC M0 = M[n-1] = 0
                M = new double[n];
                M[0] = 0.0;
                for (int i = 1; i <= n - 2; ++i) M[i] = Mint[i - 1];
                M[n - 1] = 0.0;
            }

            // compute per-interval coefficients A,B,C,D
            A = new double[m];
            B = new double[m];
            C = new double[m];
            D = new double[m];
            for (int i = 0; i < m; ++i)
            {
                double hh = h[i];
                A[i] = y[i];
                // B = (y[i+1]-y[i])/h - h*(2*M[i] + M[i+1]) / 6
                B[i] = (y[i + 1] - y[i]) / hh - (hh * (2.0 * M[i] + M[i + 1])) / 6.0;
                C[i] = M[i] / 2.0;
                D[i] = (M[i + 1] - M[i]) / (6.0 * hh);
            }

            // prefix integrals
            prefIntegral = new double[n];
            prefIntegral[0] = 0.0;
            for (int i = 0; i < m; ++i)
            {
                double H = h[i];
                double full = A[i] * H + 0.5 * B[i] * H * H + (1.0 / 3.0) * C[i] * H * H * H + 0.25 * D[i] * H * H * H * H;
                prefIntegral[i + 1] = prefIntegral[i] + full;
            }
        }

        // evaluate S(z)
        public double Evaluate(double z)
        {
            int i = Program.binsearch(x, z);
            double t = z - x[i];
            return A[i] + B[i] * t + C[i] * t * t + D[i] * t * t * t;
        }

        // evaluate S'(z)
        public double Derivative(double z)
        {
            int i = Program.binsearch(x, z);
            double t = z - x[i];
            return B[i] + 2.0 * C[i] * t + 3.0 * D[i] * t * t;
        }

        // definite integral from x0 to z
        public double Integral(double z)
        {
            int i = Program.binsearch(x, z);
            double t = z - x[i];
            double part = A[i] * t + 0.5 * B[i] * t * t + (1.0 / 3.0) * C[i] * t * t * t + 0.25 * D[i] * t * t * t * t;
            return prefIntegral[i] + part;
        }

        // debug accessors
        public double[] GetM() { return (double[])M.Clone(); }
        public double[] GetA() { return (double[])A.Clone(); }
        public double[] GetB() { return (double[])B.Clone(); }
        public double[] GetC() { return (double[])C.Clone(); }
        public double[] GetD() { return (double[])D.Clone(); }
    }



    static void Main(string[] args)
    {
        CultureInfo CI = CultureInfo.InvariantCulture;

        // Create the tabulated data: xi = 0..9, yi = cos(xi)
        int n = 10;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; ++i)
        {
            x[i] = i;
            y[i] = Math.Cos(x[i]);
        }

        // Print a small demonstration to the terminal
        Console.WriteLine("Task A:\n");
        Console.WriteLine("Table (x_i, y_i = cos(x_i)):");
        for (int i = 0; i < n; ++i) Console.WriteLine($"  {x[i],3:0}  {y[i],12:G6}");
        Console.WriteLine();

        // Show interpolation at a few sample points (integers + midpoints)
        double[] sampleZ = { 0.0, 0.5, 1.0, 2.3, 4.7, 8.9, 9.0 };
        Console.WriteLine("Sample evaluations of linterp and linterpInteg (from 0):");
        Console.WriteLine("   z       linterp(z)     linterpInteg(z)    true integral sin(z)");
        foreach (double z in sampleZ)
        {
            double li = linterp(x, y, z);
            double integ = linterpInteg(x, y, z);
            double trueInteg = Math.Sin(z) - Math.Sin(x[0]); // sin(0)=0
            Console.WriteLine($"{z,6:0.000}  {li,13:G8}  {integ,15:G8}  {trueInteg,13:G8}");
        }
        Console.WriteLine();


        using (var w = new StreamWriter("data_points.dat"))
        {
            for (int i = 0; i < n; ++i)
            {
                w.WriteLine(string.Format(CI, "{0} {1}", x[i], y[i]));
            }
        }

        using (var w = new StreamWriter("spline.dat"))
        {
            double z0 = x[0];
            double z1 = x[n - 1];
            int M = 901; // step 0.01 -> (9-0)/0.01 = 900 steps, +1 => 901 points
            for (int k = 0; k < M; ++k)
            {
                double z = z0 + (z1 - z0) * k / (M - 1);
                double li = linterp(x, y, z);
                double integ = linterpInteg(x, y, z);
                double trueInteg = Math.Sin(z);
                w.WriteLine(string.Format(CI, "{0} {1} {2} {3}", z, li, integ, trueInteg));
            }
        }

        Console.WriteLine("To see the plots of the linear interpolation and its antiderivative along with the true values for cos(x) open the file spline.png");


        Console.WriteLine("\nTask B: Quadratic spline for y = sin(x)");

        // reuse x[] from above (xi = 0..9), but create sin-values
        double[] y_sin = new double[n];
        for (int i = 0; i < n; ++i) y_sin[i] = Math.Sin(x[i]);

        // build quadratic spline
        Qspline qs = new Qspline(x, y_sin);

        // write node file
        using (var w = new StreamWriter("qs_data_points.dat", false, System.Text.Encoding.UTF8))
        {
            for (int i = 0; i < n; ++i) w.WriteLine(string.Format(CI, "{0} {1}", x[i], y_sin[i]));
        }

        // write dense sampling file: columns: z, qs(z), qs_integral(z), trueIntegral(z)
        using (var w = new StreamWriter("qspline.dat", false, System.Text.Encoding.UTF8))
        {
            double z0 = x[0], z1 = x[n - 1];
            int M = 601;
            for (int k = 0; k < M; ++k)
            {
                double z = z0 + (z1 - z0) * k / (double)(M - 1);
                double s = qs.Evaluate(z);
                double integ = qs.Integral(z);
                double trueInteg = 1.0 - Math.Cos(z); // ∫_0^z sin(t) dt = 1 - cos(z)
                w.WriteLine(string.Format(CI, "{0} {1} {2} {3}", z, s, integ, trueInteg));
            }
        }


        // print b and c for the sin nodes for quick inspection
        double[] bCoeffs = qs.GetB();
        double[] cCoeffs = qs.GetC();
        Console.WriteLine("\nqspline coefficients for sin nodes:");
        Console.Write("b:");
        for (int i = 0; i < bCoeffs.Length; ++i) Console.Write($" {bCoeffs[i]:G6}");
        Console.WriteLine();
        Console.Write("c:");
        for (int i = 0; i < cCoeffs.Length; ++i) Console.Write($" {cCoeffs[i]:G6}");
        Console.WriteLine("\nSee the quadratic spline plots in the file qspline.png.");
        


        Console.WriteLine("\nTask C: Cubic spline for y = sin(x)");

        // reuse x[] from above (xi = 0..9), build sin-values
        double[] y_sin_c = new double[n];
        for (int i = 0; i < n; ++i) y_sin_c[i] = Math.Sin(x[i]);

        // build cubic spline (natural)
        CubicSpline cs = new CubicSpline(x, y_sin_c);

        // write nodes file
        using (var w = new StreamWriter("cubic_nodes.dat"))
        {
            for (int i = 0; i < n; ++i) w.WriteLine(string.Format(CI, "{0} {1}", x[i], y_sin_c[i]));
        }

        // write dense sampling: z, S(z), integral(z), trueIntegral
        using (var w = new StreamWriter("cubic.dat"))
        {
            double z0 = x[0], z1 = x[n - 1];
            int M = 601;
            for (int k = 0; k < M; ++k)
            {
                double z = z0 + (z1 - z0) * k / (double)(M - 1);
                double s = cs.Evaluate(z);
                double integ = cs.Integral(z);
                double trueInteg = 1.0 - Math.Cos(z); // ∫_0^z sin(t) dt
                w.WriteLine(string.Format(CI, "{0} {1} {2} {3}", z, s, integ, trueInteg));
            }
        }


        // print a few sample comparisons
        double[] sampleZc = { 0.0, 0.5, 2.0, 4.3, 8.9, 9.0 };
        Console.WriteLine("   z      cubicS(z)    cubicIntegral(z)   trueIntegral");
        foreach (double z in sampleZc)
        {
            double sv = cs.Evaluate(z);
            double iv = cs.Integral(z);
            double tv = 1.0 - Math.Cos(z);
            Console.WriteLine($"{z,6:0.000}  {sv,12:G8}  {iv,15:G8}  {tv,15:G8}");
        }

        Console.WriteLine("In the file cubic1.png we see that we can't distinguish between the cubic integral and the true integral for sin and that our cubic spline graph and the graph for the built-in cspline from gnuplot are also indistiguishable, so it seems that gnuplot indeed produces a similar cubic spline to our implementation.");



        // End
    }
}
