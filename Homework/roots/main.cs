// main.cs
// Newton's method (multivariate) with numerical Jacobian and backtracking line-search.
// Uses the QR class (modified Gram-Schmidt) you supplied.
// Compile with: mcs main.cs -out:main.exe   (or use your makefile / dotnet as appropriate)

using System;
using System.Text;
using System.IO;
using System.Globalization;
using System.Collections.Generic;

namespace QRDecompositionExample
{
    // --------------------------
    // Simple vector class
    // --------------------------
    public class vector
    {
        public double[] data;
        public int size => data.Length;

        public double this[int i]
        {
            get => data[i];
            set => data[i] = value;
        }

        public vector(int n) { data = new double[n]; }

        public vector(double[] arr) { data = (double[])arr.Clone(); }

        public vector copy() { return new vector((double[])data.Clone()); }

        public override string ToString()
        {
            return "(" + string.Join(", ", Array.ConvertAll(data, d => d.ToString("G6", CultureInfo.InvariantCulture))) + ")";
        }

        public double norm()
        {
            double s = 0.0;
            for (int i = 0; i < size; i++) s += data[i] * data[i];
            return Math.Sqrt(s);
        }

        public static vector operator +(vector a, vector b)
        {
            if (a.size != b.size) throw new Exception("vector size mismatch");
            var c = new vector(a.size);
            for (int i = 0; i < a.size; i++) c[i] = a[i] + b[i];
            return c;
        }

        public static vector operator -(vector a, vector b)
        {
            if (a.size != b.size) throw new Exception("vector size mismatch");
            var c = new vector(a.size);
            for (int i = 0; i < a.size; i++) c[i] = a[i] - b[i];
            return c;
        }

        public static vector operator *(double s, vector a)
        {
            var c = new vector(a.size);
            for (int i = 0; i < a.size; i++) c[i] = s * a[i];
            return c;
        }

        public static vector operator *(vector a, double s) => s * a;

        public static vector operator -(vector a) {
            var c = new vector(a.size);
            for (int i = 0; i < a.size; i++) c[i] = -a[i];
            return c;
        }
    }

    // --------------------------
    // Simple matrix class (column-major)
    // --------------------------
    public class matrix
    {
        public readonly int size1, size2; // rows, cols
        private double[] data; // column-major: data[i + j*size1]

        public matrix(int n, int m)
        {
            size1 = n; size2 = m;
            data = new double[n * m];
        }

        // convenience constructor for square matrix
        public matrix(int n) : this(n, n) { }

        public double this[int i, int j]
        {
            get => data[i + j * size1];
            set => data[i + j * size1] = value;
        }

        public matrix copy()
        {
            matrix M = new matrix(size1, size2);
            Array.Copy(data, M.data, data.Length);
            return M;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            for (int i = 0; i < size1; i++)
            {
                for (int j = 0; j < size2; j++)
                {
                    sb.Append(this[i, j].ToString("F6", CultureInfo.InvariantCulture));
                    if (j < size2 - 1) sb.Append("\t");
                }
                sb.AppendLine();
            }
            return sb.ToString();
        }

        public static matrix Transpose(matrix M)
        {
            matrix T = new matrix(M.size2, M.size1);
            for (int i = 0; i < M.size1; i++)
                for (int j = 0; j < M.size2; j++)
                    T[j, i] = M[i, j];
            return T;
        }
    }

    // --------------------------
    // QR decomposition (modified Gram-Schmidt) from your earlier script
    // --------------------------
    public class QR
    {
        public matrix Q; // n x m
        public matrix R; // m x m

        // Constructor: perform QR on A (n x m)
        public QR(matrix A)
        {
            int n = A.size1;
            int m = A.size2;
            Q = A.copy();
            R = new matrix(m, m);

            for (int j = 0; j < m; j++)
            {
                double norm = 0.0;
                for (int i = 0; i < n; i++) norm += Q[i, j] * Q[i, j];
                norm = Math.Sqrt(norm);
                R[j, j] = norm;

                if (norm == 0.0) throw new Exception("Zero column encountered in QR (matrix rank deficient?)");

                for (int i = 0; i < n; i++) Q[i, j] /= norm;

                for (int k = j + 1; k < m; k++)
                {
                    double dot = 0.0;
                    for (int i = 0; i < n; i++) dot += Q[i, j] * Q[i, k];
                    R[j, k] = dot;
                    for (int i = 0; i < n; i++) Q[i, k] -= Q[i, j] * dot;
                }
            }
        }
        public vector solve(vector b)
        {
            int n = Q.size1;
            int m = Q.size2;
            if (b.size != n) throw new Exception("Dimension mismatch in QR.solve");

            vector y = new vector(m);
            for (int j = 0; j < m; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++) sum += Q[i, j] * b[i];
                y[j] = sum;
            }

            // back substitution R x = y
            vector x = new vector(m);
            for (int j = m - 1; j >= 0; j--)
            {
                double s = y[j];
                for (int k = j + 1; k < m; k++) s -= R[j, k] * x[k];
                x[j] = s / R[j, j];
            }
            return x;
        }

        // multiply matrices (utility)
        public static matrix Multiply(matrix A, matrix B)
        {
            if (A.size2 != B.size1) throw new Exception("Incompatible dimensions for multiplication.");
            matrix C = new matrix(A.size1, B.size2);
            for (int i = 0; i < A.size1; i++)
                for (int j = 0; j < B.size2; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < A.size2; k++) sum += A[i, k] * B[k, j];
                    C[i, j] = sum;
                }
            return C;
        }

        // determinant of R (upper triangular)
        public double det()
        {
            int m = R.size1;
            double p = 1.0;
            for (int i = 0; i < m; i++) p *= R[i, i];
            return p;
        }
    }

    class Program
    {
        // -------------------- ODE helpers (use the existing 'vector' class) --------------------
// Embedded pair: Euler(1) / Midpoint(2) like in your single-file driver
public static (vector, vector) rkstep12_ode(Func<double, vector, vector> f, double x, vector y, double h)
{
    vector k0 = f(x, y);
    vector k1 = f(x + h / 2.0, y + k0 * (h / 2.0));
    vector yh = y + (k1 * h);
    vector dy = (k1 - k0) * h;
    return (yh, dy);
}

// Adaptive driver returns lists of accepted x and y (similar to your driver)
public static (List<double>, List<vector>) driver_ode(
    Func<double, vector, vector> F,
    (double, double) interval,
    vector yinit,
    double h = 0.125,
    double acc = 1e-6,
    double eps = 1e-6
)
{
    var (a, b) = interval;
    double x = a;
    vector y = yinit.copy();
    var xlist = new List<double>() { x };
    var ylist = new List<vector>() { y.copy() };

    double intervalLength = b - a;
    if (intervalLength <= 0) return (xlist, ylist);

    double safety = 0.95;
    double maxIncrease = 2.0;
    double minStep = 1e-14;

    int p = 2; // midpoint is order 2
    double exponent = 1.0 / (p + 1.0); // 1/3

    while (true)
    {
        if (x >= b - 1e-14) return (xlist, ylist);
        if (x + h > b) h = b - x;
        if (h < minStep) throw new Exception("Step size underflow in driver_ode.");

        var (yh, dY) = rkstep12_ode(F, x, y, h);
        double tol = (acc + eps * yh.norm()) * Math.Sqrt(h / intervalLength);
        double err = dY.norm();

        if (err <= tol)
        {
            x += h;
            y = yh;
            xlist.Add(x);
            ylist.Add(y.copy());
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

// Integrate on a given r-grid and return values exactly at the r-points
public static (double[], vector[]) integrate_on_grid_ode(
    Func<double, vector, vector> F,
    double[] times,
    vector y0,
    double initialH = 0.125,
    double acc = 1e-6,
    double eps = 1e-6
)
{
    int n = times.Length;
    var ys = new vector[n];
    double currentT = times[0];
    vector currentY = y0.copy();
    ys[0] = currentY.copy();

    for (int i = 1; i < n; i++)
    {
        double tnext = times[i];
        var (xs, yslist) = driver_ode(F, (currentT, tnext), currentY, initialH, acc, eps);
        currentY = yslist[yslist.Count - 1].copy();
        currentT = tnext;
        ys[i] = currentY.copy();
    }

    return (times, ys);
}

// small Linspace convenience
public static double[] LinspaceDouble(double a, double b, int n)
{
    var t = new double[n];
    for (int i = 0; i < n; i++) t[i] = a + (b - a) * i / (n - 1);
    return t;
}

        static double rosen_scalar(vector v)
{
    double x = v[0], y = v[1];
    return (1.0 - x)*(1.0 - x) + 100.0 * (y - x*x)*(y - x*x);
}

static double himmelblau_scalar(vector v)
{
    double x = v[0], y = v[1];
    double a = x*x + y - 11.0;
    double b = x + y*y - 7.0;
    return a*a + b*b;
}
        public static matrix jacobian(Func<vector, vector> f, vector x, vector fx = null, vector dx = null)
        {
            int n = x.size;
            if (dx == null)
            {
                dx = new vector(n);
                double eps = Math.Pow(2.0, -26); 
                for (int i = 0; i < n; i++)
                    dx[i] = Math.Max(Math.Abs(x[i]), 1.0) * eps;
            }
            if (fx == null) fx = f(x);

            matrix J = new matrix(n, n); 
            for (int j = 0; j < n; j++)
            {
                x[j] += dx[j];
                vector fpx = f(x);
                for (int i = 0; i < n; i++)
                    J[i, j] = (fpx[i] - fx[i]) / dx[j];
                x[j] -= dx[j];
            }
            return J;
        }

        public static vector newton(Func<vector, vector> f, vector start, double acc = 1e-2, vector dx = null, int maxIter = 200)
        {
            double lambdaMin = 1e-12;

            vector x = start.copy();
            int n = x.size;
            if (dx == null)
            {
                dx = new vector(n);
                double eps = Math.Pow(2.0, -26);
                for (int i = 0; i < n; i++) dx[i] = Math.Max(Math.Abs(x[i]), 1.0) * eps;
            }

            vector fx = f(x);
            int iter = 0;
            while (true)
            {
                double fnorm = fx.norm();
                if (fnorm < acc) // success
                {
                    Console.WriteLine($"[newton] converged: ||f|| = {fnorm:E6} in {iter} iterations");
                    return x;
                }
                if (iter >= maxIter)
                {
                    Console.WriteLine($"[newton] reached maxIter={maxIter}, ||f|| = {fnorm:E6}");
                    return x;
                }

                // Build Jacobian at x
                matrix J = jacobian(f, x, fx, dx);

                // Solve J * Dx = -fx  (using QR)
                vector negFx = -fx;
                QR qr;
                vector Dx;
                try
                {
                    qr = new QR(J);
                    Dx = qr.solve(negFx);
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"[newton] linear solve failed at iter {iter}: {ex.Message}");
                    return x;
                }

                // If the step is already smaller than dx (componentwise), we stop: no useful progress
                bool smallStep = true;
                for (int j = 0; j < n; j++)
                {
                    if (Math.Abs(Dx[j]) >= Math.Abs(dx[j]) )
                    {
                        smallStep = false; break;
                    }
                }
                if (smallStep)
                {
                    Console.WriteLine($"[newton] step smaller than dx at iter {iter}. ||f|| = {fnorm:E6}");
                    return x;
                }

                // Backtracking line search: try lambda=1 then halve until condition satisfied
                double lambda = 1.0;
                vector z = null;
                vector fz = null;
                while (true)
                {
                    // z = x + lambda * Dx
                    z = x + (lambda * Dx);
                    fz = f(z);
                    double fznorm = fz.norm();

                    // Armijo-like condition suggested in assignment:
                    // fz.norm() < (1 - lambda/2) * fx.norm()
                    if (fznorm < (1.0 - lambda / 2.0) * fnorm)
                    {
                        break; // accept
                    }
                    if (lambda < lambdaMin)
                    {
                        // cannot decrease lambda more: accept last computed z (or keep original x)
                        // we'll accept z anyway to avoid infinite loop, but warn
                        // but if fz is worse than fx we might prefer to not update; still follow instruction.
                        // Accept z.
                        Console.WriteLine($"[newton] warning: lambda < lambdaMin at iter {iter}. Accepting current z.");
                        break;
                    }
                    lambda /= 2.0;
                }

                // After line search, check if actual step is smaller than dx (componentwise)
                bool stepTooSmall = true;
                for (int j = 0; j < n; j++)
                {
                    if (Math.Abs(lambda * Dx[j]) >= Math.Abs(dx[j]))
                    {
                        stepTooSmall = false; break;
                    }
                }
                if (stepTooSmall)
                {
                    Console.WriteLine($"[newton] step after line-search smaller than dx at iter {iter}. ||f|| = {fnorm:E6}");
                    return x;
                }

                // update
                x = z;
                fx = fz;
                iter++;

                // Optional: small info each major iteration (comment out if verbose)
                Console.WriteLine($"iter {iter,3}: ||f|| = {fx.norm():E6}, lambda = {lambda:E3}, ||Dx|| = {Dx.norm():E6}");
            }
        }

        // --------------------------
        // Test functions and gradients
        // --------------------------

        // 1D test: f(x) = x^2 - 2  (root sqrt(2))
        static vector f1d(vector v)
        {
            vector r = new vector(1);
            double x = v[0];
            r[0] = x * x - 2.0;
            return r;
        }

        // Rosenbrock function gradient (analytic). F(x,y) = (1-x)^2 + 100 (y - x^2)^2
        // grad = [ -2(1-x) - 400 x (y - x^2) , 200 (y - x^2) ]
        static vector grad_rosen(vector v)
        {
            var g = new vector(2);
            double x = v[0], y = v[1];
            g[0] = -2.0 * (1.0 - x) - 400.0 * x * (y - x * x);
            g[1] = 200.0 * (y - x * x);
            return g;
        }

        // Himmelblau gradient (analytic).
        // G(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
        // grad_x = 4x (x^2 + y - 11) + 2 (x + y^2 - 7)
        // grad_y = 2 (x^2 + y - 11) + 4y (x + y^2 - 7)
        static vector grad_himmelblau(vector v)
        {
            var g = new vector(2);
            double x = v[0], y = v[1];
            double a = x * x + y - 11.0;
            double b = x + y * y - 7.0;
            g[0] = 4.0 * x * a + 2.0 * b;
            g[1] = 2.0 * a + 4.0 * y * b;
            return g;
        }
        // put this near your other helper functions:
static vector f2d(vector v) {
    var r = new vector(2);
    r[0] = v[0]*v[0] + v[1] - 2;
    r[1] = v[0] + v[1]*v[1] - 2;
    return r;
}

        static void PrintResult(string label, vector x, Func<vector, vector> f, Func<vector,double> scalarFunc = null)
{
    Console.WriteLine($"--- {label} ---");
    Console.WriteLine($"x = {x}");
    var fx = f(x);
    Console.WriteLine($"||f(x)|| = {fx.norm():E6}");  
    if (scalarFunc != null)
        Console.WriteLine($"scalar f(x) = {scalarFunc(x):E12}");
    Console.WriteLine();
}

        static void Main(string[] args)
        {
            Console.WriteLine("Task A:\nWe start by debugging our root-finding routine using the simple one- and two-dimensional equations f(x) = x^2 - 2 and 2D test: f(x,y) = [x^2 + y - 2,  x + y^2 - 2");

            // 1) Simple 1D test: root of x^2 - 2
            {
                Console.WriteLine("1D test: f(x) = x^2 - 2");
                var start = new vector(1); start[0] = 1.0;
                var dx = new vector(1); dx[0] = Math.Max(Math.Abs(start[0]), 1.0) * Math.Pow(2.0, -26);
                vector root = newton(f1d, start, 1e-10, dx, maxIter: 50);
                PrintResult("1D result", root, f1d);
            }
// 2) Simple nonlinear 2D test: root at (1,1)
{
    Console.WriteLine("2D test: f(x,y) = [x^2 + y - 2,  x + y^2 - 2]");
    vector start = new vector(new double[]{2.0, 2.0});
    vector dx = new vector(2);
    double eps = Math.Pow(2.0, -26);
    for (int i = 0; i < 2; i++) dx[i] = Math.Max(Math.Abs(start[i]), 1.0) * eps;
    vector root = newton(f2d, start, 1e-10, dx, maxIter: 100);
    PrintResult("2D nonlinear result", root, f2d);
}

            Console.WriteLine("\nNow trying to find the extremums of the Rosenbrock valley function, which according to the linked wikipedia article is only at (1,1):\n");
            {
                Console.WriteLine("Rosenbrock: find stationary point of F(x,y) by solving grad F = 0");
                vector start = new vector(2); start[0] = 1.1; start[1] = 1.1; // classic start
                vector dx = new vector(2);
                double eps = Math.Pow(2.0, -26);
                for (int i = 0; i < 2; i++) dx[i] = Math.Max(Math.Abs(start[i]), 1.0) * eps;
                vector root = newton(grad_rosen, start, 1e-8, dx, maxIter: 200);
                PrintResult("Rosenbrock result", root, grad_rosen, rosen_scalar);
            }

            // 3) Himmelblau: multiple minima. Try several starts.
            {
                Console.WriteLine("Himmelblau: find minima by solving grad = 0 from different starts. Here there should be four local minima and one local maximum according to Wiki:");
                var starts = new vector[]
                {
                    new vector(new double[]{ 3.0, 2.0 }),
                    new vector(new double[]{ -2.0, 3.0 }),
                    new vector(new double[]{ -3.0, -3.0 }),
                    new vector(new double[]{ 3.0, -2.0 }),
                    new vector(new double[]{ 0.0, 0.0 })
                };
                int k = 0;
                foreach (var s in starts)
                {
                    k++;
                    Console.WriteLine($"\nHimmelblau start #{k}: {s}");
                    vector dx = new vector(2);
                    double eps = Math.Pow(2.0, -26);
                    for (int i = 0; i < 2; i++) dx[i] = Math.Max(Math.Abs(s[i]), 1.0) * eps;
                    vector root = newton(grad_himmelblau, s, 1e-8, dx, maxIter: 200);
                    PrintResult($"Himmelblau result #{k}", root, grad_himmelblau, himmelblau_scalar);
                }

            }

            Console.WriteLine("We have found all the local minima that match the wiki article linked (It said minima at approximately (3,2), (-2.8,3.1), (-3.8, -3.3), (3.6, -1.8) and maxima at (-0.3,-0.9))\nWe can tease out whether an extremum was min or max from the 'scalar f(x)' value in this case.\nFor both Rosenbrock's valley function and Himmelblau's function we found the extrema by first calculating the partial derivatives analytically and then searching for their roots.");

            // -------------------- Task B: Hydrogen s-wave shooting method --------------------
Console.WriteLine();
Console.WriteLine("Task B: Hydrogen s-wave shooting (shooting + bisection).");

// Parameters you can tweak:
double rmin_default = 1e-6;
double rmax_default = 8.0;
double initialH = 0.01;
double accTol = 1e-6;
double relTol = 1e-6;

// ODE: y[0]=f, y[1]=f'
Func<double, vector, double, vector> radialSystemFactory = (double E, vector dummy, double unused) =>
{
    // placeholder to allow the closure below — we'll actually make a Func in shoot_M
    return null;
};

// Shooting function M(E): integrate radial eqn and return f(rmax)
Func<double, double, double, double, double, double> shoot_M = (double E, double rmin, double rmax, double acc, double eps) =>
{
    // Initial conditions from series: f(r) ~ r - r^2 , f'(r) ~ 1 - 2r
    vector y0 = new vector(2);
    y0[0] = rmin - rmin * rmin;
    y0[1] = 1.0 - 2.0 * rmin;

    Func<double, vector, vector> F = (r, y) =>
    {
        vector dy = new vector(2);
        dy[0] = y[1];
        // avoid division by zero; r >= rmin > 0
        dy[1] = -2.0 * (E + 1.0 / r) * y[0];
        return dy;
    };

    // integrate from rmin to rmax using driver_ode (we only need final y)
    try
    {
        var (xs, ys) = driver_ode(F, (rmin, rmax), y0, initialH, acc, eps);
        vector yfinal = ys[ys.Count - 1];
        return yfinal[0];
    }
    catch (Exception ex)
    {
        // If integration failed (step underflow, blowup) return a large value with sign
        Console.WriteLine($"[shoot_M] integration failed for E={E:E}: {ex.Message}");
        return double.PositiveInfinity;
    }
};

// Helper: try to find a bracket [Ea, Eb] where M(Ea)*M(Eb) < 0 by scanning
Func<Func<double,double>, double, double, int, (bool, double, double)> find_bracket =
    (Func<double,double> Mfun, double Emin, double Emax, int Nsteps) =>
{
    double a = Emin;
    double Fa = Mfun(a);
    double step = (Emax - Emin) / Nsteps;
    for (int i = 1; i <= Nsteps; i++)
    {
        double b = Emin + i * step;
        double Fb = Mfun(b);
        if (double.IsInfinity(Fa) || double.IsInfinity(Fb))
        {
            a = b;
            Fa = Fb;
            continue;
        }
        if (Fa * Fb <= 0.0)
        {
            return (true, a, b);
        }
        a = b; Fa = Fb;
    }
    return (false, 0.0, 0.0);
};

// Bisection solver for M(E)=0 on a bracket
Func<Func<double,double>, double, double, double, int, double> bisection =
    (Func<double,double> Mfun, double a, double b, double tolE, int maxIter) =>
{
    double Fa = Mfun(a);
    double Fb = Mfun(b);
    if (double.IsInfinity(Fa) || double.IsInfinity(Fb)) throw new Exception("Infinite M at bracket endpoints.");
    if (Fa * Fb > 0.0) throw new Exception("No sign change in initial bracket.");

    double mid = 0;
    for (int iter = 0; iter < maxIter; iter++)
    {
        mid = 0.5 * (a + b);
        double Fm = Mfun(mid);
        if (double.IsInfinity(Fm)) // treat as failure; narrow bracket slightly
        {
            a = mid; continue;
        }
        if (Math.Abs(Fm) < 1e-12 || (b - a) / 2.0 < tolE) return mid;
        if (Fa * Fm <= 0.0)
        {
            b = mid; Fb = Fm;
        }
        else
        {
            a = mid; Fa = Fm;
        }
    }
    return mid;
};

// -------------- Run: find ground-state energy ----------------
double E_min = -1.5, E_max = -0.01;
int scanN = 50;

Console.WriteLine("Scanning for brackets in E in [{0}, {1}] ...", E_min, E_max);

// Wrap shoot_M into a single-argument function with default parameters
Func<double,double> Mwrapped = (double E) => shoot_M(E, rmin_default, rmax_default, accTol, relTol);

var (found, Ea, Eb) = find_bracket(Mwrapped, E_min, E_max, scanN);
if (!found)
{
    Console.WriteLine("No sign change found in the scan. Try wider E-range or check integration settings.");
}
else
{
    Console.WriteLine($"Bracket found: [{Ea:E12}, {Eb:E12}]");
    double tolE = 1e-10;
    double Eroot = bisection(Mwrapped, Ea, Eb, tolE, 100);
    Console.WriteLine($"Found root E0 = {Eroot:E12} (exact -0.5)");

    // Integrate once with Eroot on a dense r-grid and write wavefunction to file
    int Nr = 1000;
    var rgrid = LinspaceDouble(rmin_default, rmax_default, Nr);
    vector y0_final = new vector(2);
    y0_final[0] = rmin_default - rmin_default * rmin_default;
    y0_final[1] = 1.0 - 2.0 * rmin_default;
    Func<double, vector, vector> Ffinal = (r, y) =>
    {
        vector dy = new vector(2);
        dy[0] = y[1];
        dy[1] = -2.0 * (Eroot + 1.0 / r) * y[0];
        return dy;
    };
    var (rs_out, ys_out) = integrate_on_grid_ode(Ffinal, rgrid, y0_final, initialH: initialH, acc: accTol, eps: relTol);

    // Normalize the wavefunction (L2 on r-grid) if desired before writing
    // simple trapezoidal integral for norm of f(r): I = integral f(r)^2 dr
    double I = 0.0;
    for (int i = 0; i < rs_out.Length - 1; i++)
    {
        double dr = rs_out[i + 1] - rs_out[i];
        double f1 = ys_out[i][0];
        double f2 = ys_out[i + 1][0];
        I += 0.5 * (f1 * f1 + f2 * f2) * dr;
    }
    double normFactor = Math.Sqrt(I);
    if (normFactor == 0.0) normFactor = 1.0;

    string wavefile = "wavefunc.dat";
    using (var w = new StreamWriter(wavefile))
    {
        for (int i = 0; i < rs_out.Length; i++)
        {
            double r = rs_out[i];
            double f = ys_out[i][0] / normFactor; // normalized reduced radial function
            double analytic = r * Math.Exp(-r);   // shape to compare (unnormalized)
            w.WriteLine($"{r:E12} {f:E12} {analytic:E12}");
        }
    }
    Console.WriteLine($"Wavefunction written to {wavefile}. \nSee plot of the wavefunction in the file wavefunc.png");

    // ----------------- Convergence studies -> write dat files -----------------
    // 1) Vary rmin (keep others fixed)
    double[] rminList = new double[] { 1e-6, 1e-5, 1e-4, 1e-3 };
    using (var w = new StreamWriter("conv_rmin.dat"))
    {
        w.WriteLine("# rmin   E0   error(E0+0.5)");
        foreach (var rmin in rminList)
        {
            Func<double,double> Mwrap2 = (double E) => shoot_M(E, rmin, rmax_default, accTol, relTol);
            var (ok, a2, b2) = find_bracket(Mwrap2, E_min, E_max, scanN);
            if (!ok) { w.WriteLine($"{rmin:E12} NaN NaN"); continue; }
            double Eroot2 = bisection(Mwrap2, a2, b2, 1e-10, 100);
            double err = Math.Abs(Eroot2 + 0.5);
            w.WriteLine($"{rmin:E12} {Eroot2:E12} {err:E12}");
        }
    }
    Console.WriteLine("Convergence data (varying rmin) written to conv_rmin.dat\nSee the plot in the file conv_rmax.png");


    // 2) Vary rmax
    double[] rmaxList = new double[] { 4.0, 6.0, 8.0, 12.0, 16.0 };
    using (var w = new StreamWriter("conv_rmax.dat"))
    {
        w.WriteLine("# rmax   E0   error(E0+0.5)");
        foreach (var rmax in rmaxList)
        {
            Func<double,double> Mwrap3 = (double E) => shoot_M(E, rmin_default, rmax, accTol, relTol);
            var (ok, a3, b3) = find_bracket(Mwrap3, E_min, E_max, scanN);
            if (!ok) { w.WriteLine($"{rmax:E12} NaN NaN"); continue; }
            double Eroot3 = bisection(Mwrap3, a3, b3, 1e-10, 100);
            double err = Math.Abs(Eroot3 + 0.5);
            w.WriteLine($"{rmax:E12} {Eroot3:E12} {err:E12}");
        }
    }
    Console.WriteLine("Convergence data (varying rmax) written to conv_rmax.dat\nPlot can be seen in the file conv_rmax.png");

    // 3) Vary tolerances (acc/eps)
    double[] tolList = new double[] { 1e-3, 1e-4, 1e-5, 1e-6 };
    using (var w = new StreamWriter("conv_tol.dat"))
    {
        w.WriteLine("# tol   E0   error(E0+0.5)");
        foreach (var tol in tolList)
        {
            Func<double,double> Mwrap4 = (double E) => shoot_M(E, rmin_default, rmax_default, tol, tol);
            var (ok, a4, b4) = find_bracket(Mwrap4, E_min, E_max, scanN);
            if (!ok) { w.WriteLine($"{tol:E12} NaN NaN"); continue; }
            double Eroot4 = bisection(Mwrap4, a4, b4, 1e-10, 100);
            double err = Math.Abs(Eroot4 + 0.5);
            w.WriteLine($"{tol:E12} {Eroot4:E12} {err:E12}");
        }
    }
    Console.WriteLine("Convergence data (varying tolerances) written to conv_tol.dat.\nPlot in conv_tol.png.");
}

        }
    }
}

