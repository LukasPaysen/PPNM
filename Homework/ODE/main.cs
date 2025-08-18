
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

// Single-file RK solver: rkstep12 + adaptive driver + three example problems.
// Writes data file for the first example (u'' = -u).

public class Vector
{
    public double[] v;
    public int Length => v.Length;
    public Vector(int n) { v = new double[n]; }
    public Vector(params double[] arr) { v = (double[])arr.Clone(); }
    public Vector Copy() => new Vector((double[])v.Clone());
    public static Vector operator +(Vector a, Vector b)
    {
        var r = new Vector(a.Length);
        for (int i = 0; i < a.Length; i++) r.v[i] = a.v[i] + b.v[i];
        return r;
    }
    public static Vector operator -(Vector a, Vector b)
    {
        var r = new Vector(a.Length);
        for (int i = 0; i < a.Length; i++) r.v[i] = a.v[i] - b.v[i];
        return r;
    }
    public static Vector operator *(Vector a, double s)
    {
        var r = new Vector(a.Length);
        for (int i = 0; i < a.Length; i++) r.v[i] = a.v[i] * s;
        return r;
    }
    public double Norm()
    {
        double s = 0;
        for (int i = 0; i < Length; i++) s += v[i] * v[i];
        return Math.Sqrt(s);
    }
    public override string ToString() => "[" + string.Join(", ", v.Select(x => x.ToString("G6"))) + "]";
    public static Vector FromArray(double[] arr) => new Vector((double[])arr.Clone());
}

class Program
{
    // rkstep12: embedded Euler (order1) / Midpoint (order2) pair
    public static (Vector, Vector) rkstep12(Func<double, Vector, Vector> f, double x, Vector y, double h)
    {
        Vector k0 = f(x, y);
        Vector k1 = f(x + h / 2.0, y + k0 * (h / 2.0));
        Vector yh = y + k1 * h;
        Vector dy = (k1 - k0) * h;
        return (yh, dy);
    }

    // Adaptive driver (returns lists of accepted x and y)
    public static (List<double>, List<Vector>) driver(
        Func<double, Vector, Vector> F,
        (double, double) interval,
        Vector yinit,
        double h = 0.125,
        double acc = 0.01,
        double eps = 0.01
    )
    {
        var (a, b) = interval;
        double x = a;
        Vector y = yinit.Copy();
        var xlist = new List<double>();
        var ylist = new List<Vector>();
        xlist.Add(x);
        ylist.Add(y.Copy());

        double intervalLength = b - a;
        if (intervalLength <= 0) return (xlist, ylist);

        double safety = 0.95;
        double maxIncrease = 2.0;
        double minStep = 1e-14;

        int p = 2; // midpoint has order 2
        double exponent = 1.0 / (p + 1.0); // 1/(p+1) => 1/3

        while (true)
        {
            if (x >= b - 1e-14) return (xlist, ylist);
            if (x + h > b) h = b - x;
            if (h < minStep) throw new Exception("Step size underflow.");

            var (yh, dY) = rkstep12(F, x, y, h);
            double tol = (acc + eps * yh.Norm()) * Math.Sqrt(h / intervalLength);
            double err = dY.Norm();

            if (err <= tol)
            {
                x += h;
                y = yh;
                xlist.Add(x);
                ylist.Add(y.Copy());
            }

            if (err > 0.0)
            {
                double factor = Math.Min(Math.Pow(tol / err, exponent) * safety, maxIncrease);
                double minDecrease = 0.1;
                factor = Math.Max(factor, minDecrease);
                h *= factor;
            }
            else
            {
                h *= 2.0;
            }
        }
    }

    // Integrate on a given time grid, returning values exactly at the times[]
    public static (double[], Vector[]) integrate_on_grid(
        Func<double, Vector, Vector> F,
        double[] times,
        Vector y0,
        double initialH = 0.125,
        double acc = 1e-3,
        double eps = 1e-3
    )
    {
        int n = times.Length;
        var ys = new Vector[n];
        double currentT = times[0];
        Vector currentY = y0.Copy();
        ys[0] = currentY.Copy();

        for (int i = 1; i < n; i++)
        {
            double tnext = times[i];
            var (xs, yslist) = driver(F, (currentT, tnext), currentY, initialH, acc, eps);
            currentY = yslist[yslist.Count - 1].Copy();
            currentT = tnext;
            ys[i] = currentY.Copy();
        }

        return (times, ys);
    }

    public static double[] Linspace(double a, double b, int n)
    {
        var t = new double[n];
        for (int i = 0; i < n; i++) t[i] = a + (b - a) * i / (n - 1);
        return t;
    }

    static void PrintSolutionSample(string header, double[] t, Vector[] ys, int maxLines = 20)
    {
        Console.WriteLine(header);
        int n = t.Length;
        if (n <= maxLines)
        {
            for (int i = 0; i < n; i++)
                Console.WriteLine($"{t[i],7:F4} : {ys[i]}");
        }
        else
        {
            int first = 8;
            int last = 8;
            for (int i = 0; i < first; i++) Console.WriteLine($"{t[i],7:F4} : {ys[i]}");
            Console.WriteLine("   ...");
            for (int i = n - last; i < n; i++) Console.WriteLine($"{t[i],7:F4} : {ys[i]}");
        }
        Console.WriteLine();
    }

    static void Main(string[] args)
    {
        Console.WriteLine("Task A:\n");
        Console.WriteLine("We test our routine by the two ordinary differential equations u''=-u and the damped oscillator u''=-0.1u'-u and then try to recreate the example from the scipy.integrate.odeint manual. In all three cases the results are plotted");

        // -------- Example 1: u'' = -u (first-order system) --------
        Console.WriteLine();
        Console.WriteLine("Example 1: u'' = -u");
        {
            Func<double, Vector, Vector> sys1 = (t, y) => new Vector(y.v[1], -y.v[0]);
            var t1 = Linspace(0.0, 10.0, 101);
            var y0_1 = new Vector(5.0, 0.0);
            var (ts1, ys1) = integrate_on_grid(sys1, t1, y0_1, initialH: 0.05, acc: 1e-6, eps: 1e-6);
            PrintSolutionSample("Example 1 sample (t, [u, u']):", ts1, ys1);

            // Write data for gnuplot: columns -> t u u'
            string fname = "u_neg_u.dat";
            using (var w = new StreamWriter(fname))
            {
                for (int i = 0; i < ts1.Length; i++)
                {
                    double tt = ts1[i];
                    double u = ys1[i].v[0];
                    double up = ys1[i].v[1];
                    w.WriteLine($"{tt:E12} {u:E12} {up:E12}");
                }
            }
        }
        Console.WriteLine("See plot in the file u_neg_u.png");


        Console.WriteLine();
        Console.WriteLine("Example 2: damped oscillator u'' = -0.1 u' - u ");
        {
            // y[0] = u, y[1] = u'
            Func<double, Vector, Vector> sys2 = (t, y) => new Vector(y.v[1], -0.1 * y.v[1] - y.v[0]);

            var t2 = Linspace(0.0, 20.0, 201);
            var y0_2 = new Vector(1.0, 0.0);

            var (ts2, ys2) = integrate_on_grid(sys2, t2, y0_2, initialH: 0.05, acc: 1e-8, eps: 1e-8);
            PrintSolutionSample("Example 2 sample (t, [u, u']):", ts2, ys2);

            // Write data for gnuplot: columns -> t u u'
            string fname2 = "u_damped.dat";
            using (var w = new StreamWriter(fname2))
            {
                for (int i = 0; i < ts2.Length; i++)
                {
                    w.WriteLine($"{ts2[i]:E12} {ys2[i].v[0]:E12} {ys2[i].v[1]:E12}");
                }
            }
        }

        Console.WriteLine("See plot in the file u_damped.png");



        // -------- Example 3: Pendulum with friction (SciPy odeint example) --------
        Console.WriteLine();
        Console.WriteLine("Example 3: Pendulum with friction (reproducing SciPy/odeint example).");
        {
            double b = 0.25, c = 5.0;
            Func<double, Vector, Vector> pend = (t, y) =>
            {
                double theta = y.v[0], omega = y.v[1];
                return new Vector(omega, -b * omega - c * Math.Sin(theta));
            };

            var t3 = Linspace(0.0, 10.0, 101); // same grid as SciPy example
            var y0_3 = new Vector(Math.PI - 0.1, 0.0);
            var (ts3, ys3) = integrate_on_grid(pend, t3, y0_3, initialH: 0.05, acc: 1e-6, eps: 1e-6);
            PrintSolutionSample("Example 3 sample (theta, omega):", ts3, ys3);
            string fname3 = "pendulum.dat";
            using (var w = new StreamWriter(fname3))
            {
                for (int i = 0; i < ts3.Length; i++)
                {
                    w.WriteLine($"{ts3[i]:E12} {ys3[i].v[0]:E12} {ys3[i].v[1]:E12}");
                }
            }

        }

        Console.WriteLine("This plot can be seen in the file pendulum.png\n");


        Console.WriteLine("Task B:");
        // ----------------- Relativistic precession exercise -----------------
        // Add this block inside Main() where you want to run the exercise.
        // It needs Linspace(...) and integrate_on_grid(...) from your file.

        Console.WriteLine();
        Console.WriteLine("Relativistic precession of planetary orbit (u'' + u = 1 + eps u^2)");

        // common integration parameters
        double phiMax = 20.0 * Math.PI;      // several rotations (10 revolutions)
        int Nphi = 2000;                     // output sample points
        double initialH = 0.01;              // internal RK initial step (smaller for accuracy)
        double accTol = 1e-8;
        double relTol = 1e-8;


        // --- Case 1: eps=0, u(0)=1, u'(0)=0 (Newtonian circular) ---
        {
            double eps = 0.0;
            double u0 = 1.0;
            double upr0 = 0.0;
            Func<double, Vector, Vector> F = (phi, y) =>
                new Vector(y.v[1], 1.0 - y.v[0] + eps * y.v[0] * y.v[0]);

            var phiGrid = Linspace(0.0, phiMax, Nphi);
            var y0 = new Vector(u0, upr0);
            var (phis, ys) = integrate_on_grid(F, phiGrid, y0, initialH: initialH, acc: accTol, eps: relTol);

            // write phi and u to file
            string fname = "orbit_circle.dat";
            using (var w = new System.IO.StreamWriter(fname))
            {
                for (int i = 0; i < phis.Length; i++) w.WriteLine($"{phis[i]:E12} {ys[i].v[0]:E12}");
            }
            PrintSolutionSample("Relativistic exercise - Case 1 sample (phi, u):", phis, ys);
        }

        // --- Case 2: eps=0, u(0)=1, u'(0)=-0.5 (Newtonian ellipse) ---
        {
            double eps = 0.0;
            double u0 = 1.0;
            double upr0 = -0.5;
            Func<double, Vector, Vector> F = (phi, y) =>
                new Vector(y.v[1], 1.0 - y.v[0] + eps * y.v[0] * y.v[0]);

            var phiGrid = Linspace(0.0, phiMax, Nphi);
            var y0 = new Vector(u0, upr0);
            var (phis, ys) = integrate_on_grid(F, phiGrid, y0, initialH: initialH, acc: accTol, eps: relTol);

            string fname = "orbit_ellipse.dat";
            using (var w = new System.IO.StreamWriter(fname))
            {
                for (int i = 0; i < phis.Length; i++) w.WriteLine($"{phis[i]:E12} {ys[i].v[0]:E12}");
            }
            PrintSolutionSample("Relativistic exercise - Case 2 sample (phi, u):", phis, ys);
        }

        // --- Case 3: eps=0.01, u(0)=1, u'(0)=-0.5 (relativistic precession) ---
        {
            double eps = 0.01;
            double u0 = 1.0;
            double upr0 = -0.5;
            Func<double, Vector, Vector> F = (phi, y) =>
                new Vector(y.v[1], 1.0 - y.v[0] + eps * y.v[0] * y.v[0]);

            var phiGrid = Linspace(0.0, phiMax, Nphi);
            var y0 = new Vector(u0, upr0);
            var (phis, ys) = integrate_on_grid(F, phiGrid, y0, initialH: initialH, acc: accTol, eps: relTol);

            string fname = "orbit_precess.dat";
            using (var w = new System.IO.StreamWriter(fname))
            {
                for (int i = 0; i < phis.Length; i++) w.WriteLine($"{phis[i]:E12} {ys[i].v[0]:E12}");
            }
            PrintSolutionSample("Relativistic exercise - Case 3 sample (phi, u):", phis, ys);
        }

        Console.WriteLine("See the plot in pendulum.png");
        // -------------------------------------------------------------------

    }
}
