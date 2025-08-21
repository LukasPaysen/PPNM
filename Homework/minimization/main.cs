using System;
using System.Globalization;
using System.Text;

namespace QRDecompositionExample
{
    // (User-provided vector, matrix, QR classes assumed identical to what you pasted.)

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
    // QR decomposition (modified Gram-Schmidt)
    // --------------------------
    public class QR
    {
        public matrix Q; // n x m
        public matrix R; // m x m

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

            vector x = new vector(m);
            for (int j = m - 1; j >= 0; j--)
            {
                double s = y[j];
                for (int k = j + 1; k < m; k++) s -= R[j, k] * x[k];
                x[j] = s / R[j, j];
            }
            return x;
        }

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
        // --- small helper: produce a unique key for a vector state
static string VecKey(vector x) {
    var sb = new StringBuilder();
    for (int i = 0; i < x.size; ++i) {
        if (i>0) sb.Append(',');
        // Use "R" to keep a stable round-trip representation
        sb.Append(x[i].ToString("R", CultureInfo.InvariantCulture));
    }
    return sb.ToString();
}

// --- cached evaluator
static double EvalWithCache(Func<vector,double> phi, vector x, System.Collections.Generic.Dictionary<string,double> cache) {
    string key = VecKey(x);
    if (cache.TryGetValue(key, out double val)) return val;
    val = phi(x);
    cache[key] = val;
    return val;
}

// --- central-difference gradient (uses cache)
static vector GradientCentral(Func<vector,double> phi, vector x, System.Collections.Generic.Dictionary<string,double> cache) {
    int n = x.size;
    vector g = new vector(n);
    double phiX = EvalWithCache(phi, x, cache);
    for (int i = 0; i < n; ++i) {
        double dxi = (1.0 + Math.Abs(x[i])) * Math.Pow(2.0, -26);
        x[i] += dxi;
        double phip = EvalWithCache(phi, x, cache);
        x[i] -= 2.0 * dxi;
        double phim = EvalWithCache(phi, x, cache);
        x[i] += dxi;
        g[i] = (phip - phim) / (2.0 * dxi);
    }
    return g;
}

// --- central-difference Hessian via central differences of gradient (uses cache)
static matrix HessianCentral(Func<vector,double> phi, vector x) {
    int n = x.size;
    matrix H = new matrix(n, n);
    // Use a cache dictionary for phi evaluations shared across calls
    var cache = new System.Collections.Generic.Dictionary<string,double>();
    // base gradient (computed with central diffs and caching)
    vector g0 = GradientCentral(phi, x, cache);

    for (int j = 0; j < n; ++j) {
        double dxj = (1.0 + Math.Abs(x[j])) * Math.Pow(2.0, -13);
        // g(x + dx_j)
        x[j] += dxj;
        vector gplus = GradientCentral(phi, x, cache);
        // g(x - dx_j)
        x[j] -= 2.0 * dxj;
        vector gminus = GradientCentral(phi, x, cache);
        // restore
        x[j] += dxj;
        for (int i = 0; i < n; ++i) {
            H[i, j] = (gplus[i] - gminus[i]) / (2.0 * dxj);
        }
    }
    return H;
}

// --- Newton variant that uses central diffs (mirrors your Newton)
static (vector xmin, int iterations) NewtonCentral(Func<vector,double> phi, vector x0, double acc = 1e-3, int maxIter = 1000) {
    vector x = x0.copy();
    int iter = 0;
    while (iter < maxIter) {
        iter++;
        // create a fresh cache per Newton iteration (keeps memory bounded)
        var cache = new System.Collections.Generic.Dictionary<string,double>();
        double phiX = EvalWithCache(phi, x, cache);
        vector g = GradientCentral(phi, x, cache);
        if (g.norm() < acc) break;

        matrix H = HessianCentral(phi, x);
        vector negg = -g;
        vector dx = null;
        try {
            QR qr = new QR(H);
            dx = qr.solve(negg);
        } catch {
            dx = negg;
        }

        double gdotdx = 0.0;
        for (int i = 0; i < x.size; ++i) gdotdx += g[i] * dx[i];
        if (gdotdx >= 0.0) dx = negg;

        double lambda = 1.0;
        bool accepted = false;
        while (lambda >= 1.0 / 1024.0) {
            vector xtrial = new vector(x.size);
            for (int i = 0; i < x.size; ++i) xtrial[i] = x[i] + lambda * dx[i];
            double phit = EvalWithCache(phi, xtrial, cache);
            if (phit < phiX) {
                x = xtrial;
                accepted = true;
                break;
            }
            lambda *= 0.5;
        }
        if (!accepted) break;
    }
    return (x, iter);
}

        static double rosenbrock(vector v)
        {
            double x = v[0], y = v[1];
            double t = 1.0 - x;
            double s = y - x * x;
            return t * t + 100.0 * s * s;
        }

        static double himmelblau(vector v)
        {
            double x = v[0], y = v[1];
            double a = x * x + y - 11.0;
            double b = x + y * y - 7.0;
            return a * a + b * b;
        }

        static vector Gradient(Func<vector, double> phi, vector x)
        {
            int n = x.size;
            vector g = new vector(n);
            double phiX = phi(x);
            for (int i = 0; i < n; ++i)
            {
                double dxi = (1.0 + Math.Abs(x[i])) * Math.Pow(2.0, -26);
                x[i] += dxi;
                double phixdx = phi(x);
                g[i] = (phixdx - phiX) / dxi;
                x[i] -= dxi;
            }
            return g;
        }

        static matrix Hessian(Func<vector, double> phi, vector x)
        {
            int n = x.size;
            matrix H = new matrix(n, n);
            vector g0 = Gradient(phi, x);
            for (int j = 0; j < n; ++j)
            {
                double dxj = (1.0 + Math.Abs(x[j])) * Math.Pow(2.0, -13);
                x[j] += dxj;
                vector gj = Gradient(phi, x);
                for (int i = 0; i < n; ++i) H[i, j] = (gj[i] - g0[i]) / dxj;
                x[j] -= dxj;
            }
            return H;
        }

        static (vector xmin, int iterations) Newton(Func<vector, double> phi, vector x0, double acc = 1e-3, int maxIter = 1000)
        {
            vector x = x0.copy();
            int iter = 0;
            while (iter < maxIter)
            {
                iter++;
                double phiX = phi(x);
                vector g = Gradient(phi, x);
                if (g.norm() < acc) break;

                matrix H = Hessian(phi, x);
                vector negg = -g;
                vector dx = null;
                try
                {
                    QR qr = new QR(H);
                    dx = qr.solve(negg);
                }
                catch
                {
                    dx = negg; // fallback
                }

                double gdotdx = 0.0;
                for (int i = 0; i < x.size; ++i) gdotdx += g[i] * dx[i];
                if (gdotdx >= 0.0) dx = negg;

                double lambda = 1.0;
                bool accepted = false;
                while (lambda >= 1.0 / 1024.0)
                {
                    vector xtrial = new vector(x.size);
                    for (int i = 0; i < x.size; ++i) xtrial[i] = x[i] + lambda * dx[i];
                    double phit = phi(xtrial);
                    if (phit < phiX)
                    {
                        x = xtrial;
                        accepted = true;
                        break;
                    }
                    lambda *= 0.5;
                }
                if (!accepted)
                {
                    // no acceptable step found, keep current x and exit
                    break;
                }
            }
            return (x, iter);
        }

        static void Main(string[] args)
{
    Func<double,double,double,double,double> breitwigner = (E,m,G,A) => A/((E-m)*(E-m) + G*G/4.0);

    // Task A: Rosenbrock
    Console.WriteLine("\nTask A:\nFinding a minimum of Himmelblau's and Rosensbrock's valley function whilst also recording the number of steps it takes our algorithm to reach the minimum\n");
    var x0 = new vector(new double[] { -1.2, 2.0 });
    Console.WriteLine("Rosenbrock start:");
    Console.WriteLine($"x0 = {x0}");
    Console.WriteLine($"f(x0) = {rosenbrock(x0).ToString("G6", CultureInfo.InvariantCulture)}");
    var (xminR, itR) = Newton(rosenbrock, x0, 1e-6, 1000);
    Console.WriteLine("Rosenbrock result:");
    Console.WriteLine($"x = {xminR}");
    Console.WriteLine($"f(x) = {rosenbrock(xminR).ToString("G6", CultureInfo.InvariantCulture)}");
    Console.WriteLine($"iterations = {itR}");

    // Task A: Himmelblau
    var y0 = new vector(new double[] { 0.0, 0.0 });
    Console.WriteLine("\nHimmelblau start:");
    Console.WriteLine($"x0 = {y0}");
    Console.WriteLine($"f(x0) = {himmelblau(y0).ToString("G6", CultureInfo.InvariantCulture)}");
    var (xminH, itH) = Newton(himmelblau, y0, 1e-6, 1000);
    Console.WriteLine("Himmelblau result:");
    Console.WriteLine($"x = {xminH}");
    Console.WriteLine($"f(x) = {himmelblau(xminH).ToString("G6", CultureInfo.InvariantCulture)}");
    Console.WriteLine($"iterations = {itH}");

    // Task B: read data
    Console.WriteLine("\n\nTask B: \nCreating a plot of the data points along with the fit by the Breit-Wigner function.");
var energy = new System.Collections.Generic.List<double>();
var signal = new System.Collections.Generic.List<double>();
var sigmaErr = new System.Collections.Generic.List<double>();
var separators = new char[] {' ','\t'};
var options = StringSplitOptions.RemoveEmptyEntries;
string line;
while ((line = Console.In.ReadLine()) != null)
{
    if (line.Trim().Length == 0) continue;
    var words = line.Split(separators, options);
    energy.Add(double.Parse(words[0], CultureInfo.InvariantCulture));
    signal.Add(double.Parse(words[1], CultureInfo.InvariantCulture));
    sigmaErr.Add(double.Parse(words[2], CultureInfo.InvariantCulture));
}

// Deviation function
Func<vector,double> D = (vector p) => {
    double m = p[0];
    double gamma = Math.Abs(p[1]);
    double A = p[2];
    double sum = 0.0;
    for (int i = 0; i < energy.Count; ++i)
    {
        double denom = (energy[i]-m)*(energy[i]-m) + (gamma*gamma)/4.0;
        double Fi = A / denom;
        double r = (Fi - signal[i]) / sigmaErr[i];
        sum += r*r;
    }
    return sum;
};

// Fit and write outputs (use these variable names)
var start = new vector(new double[] {125, 2, 10});
var (best, itB) = Newton(D, start, 1e-6, 500);
Console.WriteLine($"Best fit: m={best[0]}, Γ={best[1]}, A={best[2]}, iterations={itB}\nSee the plot in the file higgs.png");



// Write fitted curve to fit.dat
using (var w = new System.IO.StreamWriter("fit.dat"))
{
    for (double e = 100.0; e <= 160.0; e += 0.2)
    {
        double F = breitwigner(e, best[0], Math.Abs(best[1]), best[2]);
        w.WriteLine($"{e.ToString(CultureInfo.InvariantCulture)} {F.ToString(CultureInfo.InvariantCulture)}");
    }
}

// Write original data to data.dat
using (var w = new System.IO.StreamWriter("data.dat"))
{
    for (int i = 0; i < energy.Count; i++)
    {
        w.WriteLine($"{energy[i].ToString(CultureInfo.InvariantCulture)} {signal[i].ToString(CultureInfo.InvariantCulture)} {sigmaErr[i].ToString(CultureInfo.InvariantCulture)}");
    }
}
// Compare forward vs central (do not overwrite Task A output)
Console.WriteLine("\n\nTask C:\nComparing central and forward finite difference approximation. \nWe do the exactly the same example as in task A and see if there is difference in the number of iterations.");
Console.WriteLine("\nRosenbrock start:");
Console.WriteLine($"x0 = {x0}");
Console.WriteLine($"f(x0) = {rosenbrock(x0).ToString("G6", CultureInfo.InvariantCulture)}");
Console.WriteLine("\nRosenbrock result central differences:");

var resF = Newton(rosenbrock, x0, 1e-6, 1000);
var (paa1, paa2) = NewtonCentral(rosenbrock, x0, 1e-6, 1000);

Console.WriteLine($"x = {paa1}");

Console.WriteLine($"f(x) = {rosenbrock(paa1).ToString("G6", CultureInfo.InvariantCulture)}");
Console.WriteLine($"Rosenbrock: forward iters={resF.iterations}, central iters={paa2}\nSo central differences got a slightly more precise result using two fewer iterations.");


Console.WriteLine("\nHimmelblau start:");
Console.WriteLine($"x0 = {y0}");
Console.WriteLine($"f(x0) = {himmelblau(y0).ToString("G6", CultureInfo.InvariantCulture)}");
Console.WriteLine("\nHimmelblau result central differences:");

var resEE = Newton(himmelblau, y0, 1e-6, 1000);
var (faa1, faa2) = NewtonCentral(himmelblau, y0, 1e-6, 1000);
Console.WriteLine($"x = {faa1}");

Console.WriteLine($"f(x) = {himmelblau(faa1).ToString("G6", CultureInfo.InvariantCulture)}");
Console.WriteLine($"Rosenbrock: forward iters={resEE.iterations}, central iters={faa2}\nSo central differences got a slightly more precise result but with the same amount of iterations for Himmelblau.");
/*

var startH = new vector(new double[] { 0.0, 0.0 });
var resFh = Newton(himmelblau, startH, 1e-6, 1000);
var resCh = NewtonCentral(himmelblau, startH, 1e-6, 1000);
Console.WriteLine($"Himmelblau: forward iters={resFh.iterations}, central iters={resCh.iterations}");


    
    var (xminR, itR) = Newton(rosenbrock, x0, 1e-6, 1000);
    Console.WriteLine("Rosenbrock result:");
    Console.WriteLine($"x = {xminR}");
    Console.WriteLine($"f(x) = {rosenbrock(xminR).ToString("G6", CultureInfo.InvariantCulture)}");
    Console.WriteLine($"iterations = {itR}");

    // Task A: Himmelblau
    Console.WriteLine("\nHimmelblau start:");
    Console.WriteLine($"x0 = {y0}");
    Console.WriteLine($"f(x0) = {himmelblau(y0).ToString("G6", CultureInfo.InvariantCulture)}");
    var (xminH, itH) = Newton(himmelblau, y0, 1e-6, 1000);
    Console.WriteLine("Himmelblau result:");
    Console.WriteLine($"x = {xminH}");
    Console.WriteLine($"f(x) = {himmelblau(xminH).ToString("G6", CultureInfo.InvariantCulture)}");
    Console.WriteLine($"iterations = {itH}");*/
}


    }
}
