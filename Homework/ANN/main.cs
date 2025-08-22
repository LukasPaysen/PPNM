
using System;
using System.Globalization;
using System.Text;

namespace QRDecompositionExample {

    public class vector {
        public double[] data;
        public int size => data.Length;

        public double this[int i] { get => data[i]; set => data[i] = value; }
        public vector(int n) { data = new double[n]; }
        public vector(double[] arr) { data = (double[])arr.Clone(); }
        public vector copy() { return new vector((double[])data.Clone()); }

        public override string ToString() {
            return "(" + string.Join(", ", Array.ConvertAll(data, d => d.ToString("G6", CultureInfo.InvariantCulture))) + ")";
        }

        public double norm() {
            double s = 0.0;
            for (int i = 0; i < size; i++) s += data[i] * data[i];
            return Math.Sqrt(s);
        }

        public static vector operator +(vector a, vector b) {
            if (a.size != b.size) throw new Exception("vector size mismatch");
            var c = new vector(a.size);
            for (int i = 0; i < a.size; i++) c[i] = a[i] + b[i];
            return c;
        }

        public static vector operator -(vector a, vector b) {
            if (a.size != b.size) throw new Exception("vector size mismatch");
            var c = new vector(a.size);
            for (int i = 0; i < a.size; i++) c[i] = a[i] - b[i];
            return c;
        }

        public static vector operator *(double s, vector a) {
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
    public class matrix {
        public readonly int size1, size2; // rows, cols
        private double[] data; // column-major: data[i + j*size1]
        public matrix(int n, int m) { size1 = n; size2 = m; data = new double[n * m]; }
        public matrix(int n) : this(n, n) { }

        public double this[int i, int j] { get => data[i + j * size1]; set => data[i + j * size1] = value; }

        public matrix copy() { var M = new matrix(size1, size2); Array.Copy(data, M.data, data.Length); return M; }

        public override string ToString() {
            var sb = new StringBuilder();
            for (int i = 0; i < size1; i++) {
                for (int j = 0; j < size2; j++) {
                    sb.Append(this[i, j].ToString("F6", CultureInfo.InvariantCulture));
                    if (j < size2 - 1) sb.Append("\t");
                }
                sb.AppendLine();
            }
            return sb.ToString();
        }

        public static matrix Transpose(matrix M) {
            matrix T = new matrix(M.size2, M.size1);
            for (int i = 0; i < M.size1; i++)
                for (int j = 0; j < M.size2; j++)
                    T[j, i] = M[i, j];
            return T;
        }
    }

    public class QR {
        public matrix Q; // n x m
        public matrix R; // m x m

        public QR(matrix A) {
            int n = A.size1, m = A.size2;
            Q = A.copy();
            R = new matrix(m, m);

            for (int j = 0; j < m; j++) {
                double norm = 0.0;
                for (int i = 0; i < n; i++) norm += Q[i, j] * Q[i, j];
                norm = Math.Sqrt(norm);
                if (norm == 0.0) throw new Exception("Zero column encountered in QR (matrix rank deficient?)");
                R[j, j] = norm;

                for (int i = 0; i < n; i++) Q[i, j] /= norm;

                for (int k = j + 1; k < m; k++) {
                    double dot = 0.0;
                    for (int i = 0; i < n; i++) dot += Q[i, j] * Q[i, k];
                    R[j, k] = dot;
                    for (int i = 0; i < n; i++) Q[i, k] -= Q[i, j] * dot;
                }
            }
        }
        public vector solve(vector b) {
            int n = Q.size1, m = Q.size2;
            if (b.size != n) throw new Exception("Dimension mismatch in QR.solve");
            vector y = new vector(m);
            for (int j = 0; j < m; j++) {
                double sum = 0.0;
                for (int i = 0; i < n; i++) sum += Q[i, j] * b[i];
                y[j] = sum;
            }
            vector x = new vector(m);
            for (int j = m - 1; j >= 0; j--) {
                double s = y[j];
                for (int k = j + 1; k < m; k++) s -= R[j, k] * x[k];
                x[j] = s / R[j, j];
            }
            return x;
        }
    }
        public class ann {
    public int n;
    public Func<double,double> f   = z => z * Math.Exp(-z*z);                       // wavelet
    public Func<double,double> df  = z => Math.Exp(-z*z) * (1.0 - 2.0*z*z);         // f'
    public Func<double,double> d2f = z => (4.0*z*z*z - 6.0*z) * Math.Exp(-z*z);     // f''
    public Func<double,double> fint= z => -0.5 * Math.Exp(-z*z);                    // ∫f = -½ e^{-z²}
    public vector p;

    public ann(int n){
        this.n = n;
        p = new vector(3*n); // [a | s | w]
    }
public ann(int n,
           Func<double, double> f = null,
           Func<double, double> df = null,
           Func<double, double> d2f = null,
           Func<double, double> fint = null)
{
    this.n = n;
    // f(z) = z e^{-z^2}
    this.f    = f    ?? (z => z * Math.Exp(-z * z));
    this.df   = df   ?? (z => Math.Exp(-z * z) * (1.0 - 2.0 * z * z));
    this.d2f  = d2f  ?? (z => (4.0 * z * z * z - 6.0 * z) * Math.Exp(-z * z));
    // A(z) = ∫ f(z) dz = -1/2 e^{-z^2}  (since (e^{-z^2})' = -2 z e^{-z^2})
    this.fint = fint ?? (z => -0.5 * Math.Exp(-z * z));
    p = new vector(3 * n); // [a_1..a_n | s_1..s_n | w_1..w_n],  b_i = exp(s_i)
}

public double dresponse(double x){
    if (p == null || p.size != 3*n) throw new InvalidOperationException("ANN not initialized: call train(...) first.");
    if (df == null) throw new InvalidOperationException("Activation derivative df is null.");
    double y1 = 0.0;
    for (int i = 0; i < n; i++){
        double ai = p[i], si = p[n+i], wi = p[2*n+i];
        double bi = Math.Exp(si), zi = (x - ai)/bi;
        y1 += wi * (df(zi) / bi);
    }
    return y1;
}

public double ddresponse(double x){
    if (p == null || p.size != 3*n) throw new InvalidOperationException("ANN not initialized: call train(...) first.");
    if (d2f == null) throw new InvalidOperationException("Activation second derivative d2f is null.");
    double y2 = 0.0;
    for (int i = 0; i < n; i++){
        double ai = p[i], si = p[n+i], wi = p[2*n+i];
        double bi = Math.Exp(si), zi = (x - ai)/bi;
        y2 += wi * (d2f(zi) / (bi*bi));
    }
    return y2;
}

public double primitive(double x, double c = 0.0){
    if (p == null || p.size != 3*n) throw new InvalidOperationException("ANN not initialized: call train(...) first.");
    if (fint == null) throw new InvalidOperationException("Activation antiderivative fint is null.");
    double H = 0.0;
    for (int i = 0; i < n; i++){
        double ai = p[i], si = p[n+i], wi = p[2*n+i];
        double bi = Math.Exp(si);
        double zx = (x - ai)/bi, zc = (c - ai)/bi;
        H += wi * bi * (fint(zx) - fint(zc));
    }
    return H;
}

        

        public double response(double x)
        {
            double y = 0.0;
            for (int i = 0; i < n; i++)
            {
                double ai = p[i];
                double si = p[n + i];
                double wi = p[2 * n + i];
                double bi = Math.Exp(si);
                double zi = (x - ai) / bi;
                y += wi * f(zi);
            }
            return y;
        }

        public void initialize(double[] xs, double[] ys) {
            double xmin = xs[0], xmax = xs[0];
            for (int k = 1; k < xs.Length; ++k) { if (xs[k] < xmin) xmin = xs[k]; if (xs[k] > xmax) xmax = xs[k]; }
            double L = xmax - xmin;
            if (L == 0) L = 1;

            for (int i = 0; i < n; i++) {
                double t = (n == 1) ? 0.5 : (double)i / (n - 1);
                p[i] = xmin + t * L;    // a_i
            }

            double b0 = 0.75 * (L / n);      // a sensible width
            double s0 = Math.Log(Math.Max(b0, 1e-3));
            for (int i = 0; i < n; i++) p[n + i] = s0; // s_i = ln b_i

            matrix G = new matrix(xs.Length, n);
            for (int k = 0; k < xs.Length; ++k) {
                for (int i = 0; i < n; ++i) {
                    double ai = p[i];
                    double bi = Math.Exp(p[n + i]);
                    double zi = (xs[k] - ai) / bi;
                    G[k, i] = f(zi);
                }
            }
            var qr = new QR(G);
            var yvec = new vector(ys.Length);
            for (int k = 0; k < ys.Length; ++k) yvec[k] = ys[k];
            vector w0 = qr.solve(yvec);
            for (int i = 0; i < n; i++) p[2 * n + i] = w0[i];
        }

        public (int iters, double mse) train(double[] xs, double[] ys, int maxIter = 200, double tol = 1e-8) {
            if (p == null || p.size != 3 * n) p = new vector(3 * n);
            if (xs.Length != ys.Length) throw new Exception("x/y length mismatch");

            // Initialize if first time (or keep current p if user preset)
            if (!initializedFlag) { initialize(xs, ys); initializedFlag = true; }

            int m = 3 * n;
            int N = xs.Length;
            double prevCost = double.PositiveInfinity;
            double mu = 1e-3;              // LM damping
            int iter = 0;

            // Buffers
            vector r = new vector(N);
            matrix J = new matrix(N, m);
            vector b_aug = new vector(N + m);

            while (iter < maxIter) {
                iter++;

                // Build residuals and analytic Jacobian J (N x m)
                double cost = 0.0;
                for (int k = 0; k < N; ++k) {
                    double x = xs[k];
                    double yk = ys[k];
                    double F = 0.0;

                    // Precompute per neuron
                    for (int i = 0; i < n; ++i) {
                        double ai = p[i];
                        double si = p[n + i];
                        double wi = p[2 * n + i];
                        double bi = Math.Exp(si);
                        double zi = (x - ai) / bi;
                        double fi = f(zi);
                        double dfi = df(zi);

                        F += wi * fi;

                        // dF/da_i = - (wi / b_i) * f'(z_i)
                        J[k, i] = - (wi / bi) * dfi;

                        // dF/ds_i = dF/db_i * db_i/ds_i = ( - wi * f'(z_i) * z_i / b_i ) * b_i = - wi * f'(z_i) * z_i
                        J[k, n + i] = - wi * dfi * zi;

                        // dF/dw_i = f(z_i)
                        J[k, 2 * n + i] = fi;
                    }

                    double rk = F - yk;
                    r[k] = rk;
                    cost += rk * rk;
                }
                cost *= 0.5; // standard 1/2 ||r||^2

                // Convergence checks
                if (Math.Abs(prevCost - cost) <= tol * (1 + cost)) break;
                prevCost = cost;

                // LM step: solve [J; sqrt(mu) I] * delta = -[r; 0]
                matrix A = new matrix(N + m, m);
                // Top: J
                for (int i = 0; i < N; ++i)
                    for (int j = 0; j < m; ++j)
                        A[i, j] = J[i, j];
                // Bottom: sqrt(mu) * I
                double sMu = Math.Sqrt(mu);
                for (int j = 0; j < m; ++j) A[N + j, j] = sMu;

                // RHS
                for (int i = 0; i < N; ++i) b_aug[i] = -r[i];
                for (int j = 0; j < m; ++j) b_aug[N + j] = 0.0;

                vector delta;
                try {
                    var qr = new QR(A);
                    delta = qr.solve(b_aug);
                } catch {
                    // fall back to small gradient step on weights only if J is singular
                    delta = new vector(m);
                    for (int j = 0; j < m; ++j) delta[j] = 0.0;
                    for (int i = 0; i < n; ++i) delta[2 * n + i] = -0.01 * J[0, 2 * n + i]; // tiny nudge
                }

                // Trial update with backtracking if needed
                double lambda = 1.0;
                bool accepted = false;
                vector pNew = new vector(m);
                while (lambda >= 1.0 / 10483.0) { // allow much smaller steps if needed

                for (int j = 0; j < m; ++j) pNew[j] = p[j] + lambda * delta[j];

                    // Evaluate cost at pNew (only residuals)
                    double newCost = 0.0;
                    for (int k = 0; k < N; ++k) {
                        double x = xs[k], yk = ys[k];
                        double F = 0.0;
                        for (int i = 0; i < n; ++i) {
                            double ai = pNew[i];
                            double si = pNew[n + i];
                            double wi = pNew[2 * n + i];
                            double bi = Math.Exp(si);
                            double zi = (x - ai) / bi;
                            F += wi * f(zi);
                        }
                        double rk = F - yk;
                        newCost += rk * rk;
                    }
                    newCost *= 0.5;

                    if (newCost < cost) {
                        p = pNew;
                        cost = newCost;
                        mu *= 0.3;        // relax damping on success
                        accepted = true;
                        break;
                    }
                    lambda *= 0.5;         // backtrack
                }
                if (!accepted) {
                    mu *= 7.0;             // increase damping and try next iter
                    // also stop if step is tiny
                    if (delta.norm() <= tol * (1 + p.norm())) break;
                }
            }

            // Final MSE on training set
            double mse = 0.0;
            for (int k = 0; k < xs.Length; ++k) {
                double diff = response(xs[k]) - ys[k];
                mse += diff * diff;
            }
            mse /= xs.Length;

            return (iter, mse);
        }

        private bool initializedFlag = false;
    }

    class Program {
        static string VecKey(vector x)
        {
            var sb = new StringBuilder();
            for (int i = 0; i < x.size; ++i)
            {
                if (i > 0) sb.Append(',');
                sb.Append(x[i].ToString("R", CultureInfo.InvariantCulture));
            }
            return sb.ToString();
        }

static double EvalWithCache(Func<vector,double> phi, vector x, System.Collections.Generic.Dictionary<string,double> cache) {
    string key = VecKey(x);
    if (cache.TryGetValue(key, out double val)) return val;
    val = phi(x);
    cache[key] = val;
    return val;
}

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

static matrix HessianCentral(Func<vector,double> phi, vector x) {
    int n = x.size;
    matrix H = new matrix(n, n);
    var cache = new System.Collections.Generic.Dictionary<string,double>();
    vector g0 = GradientCentral(phi, x, cache);
    for (int j = 0; j < n; ++j) {
        double dxj = (1.0 + Math.Abs(x[j])) * Math.Pow(2.0, -13);
        x[j] += dxj;
        vector gplus = GradientCentral(phi, x, cache);
        x[j] -= 2.0 * dxj;
        vector gminus = GradientCentral(phi, x, cache);
        x[j] += dxj;
        for (int i = 0; i < n; ++i) H[i, j] = (gplus[i] - gminus[i]) / (2.0 * dxj);
    }
    return H;
}

static (vector xmin, int iterations) NewtonCentral(Func<vector,double> phi, vector x0, double acc = 1e-6, int maxIter = 200) {
    vector x = x0.copy();
    int iter = 0;
    while (iter < maxIter) {
        iter++;
        var cache = new System.Collections.Generic.Dictionary<string,double>();
        double phiX = EvalWithCache(phi, x, cache);
        vector g = GradientCentral(phi, x, cache);
        if (g.norm() < acc) break;

        matrix H = HessianCentral(phi, x);
        vector negg = -g, dx;
        try { dx = new QR(H).solve(negg); } catch { dx = negg; }

        double gdotdx = 0.0;
        for (int i = 0; i < x.size; ++i) gdotdx += g[i] * dx[i];
        if (gdotdx >= 0.0) dx = negg;

        double lambda = 1.0;
        bool accepted = false;
        while (lambda >= 1.0 / 1024.0) {
            vector xtrial = new vector(x.size);
            for (int i = 0; i < x.size; ++i) xtrial[i] = x[i] + lambda * dx[i];
            double phit = EvalWithCache(phi, xtrial, cache);
            if (phit < phiX) { x = xtrial; accepted = true; break; }
            lambda *= 0.5;
        }
        if (!accepted) break;
    }
    return (x, iter);
}

        // Example Phi: y'' + y = 0
static double Phi(double y2, double y1, double y0, double x) => y2 + y0;

static (int iters, double cost) TrainODE(ann net, double a, double b, int M,
                                         double c, double yc, double ypc)
{
    double[] xg = new double[M];
    double[] wg = new double[M]; // trapezoid weights
    for (int j = 0; j < M; ++j) {
        xg[j] = a + (b-a) * j/(M-1.0);
        wg[j] = (j==0 || j==M-1) ? 0.5*(b-a)/(M-1.0) : (b-a)/(M-1.0);
    }

    Func<vector,double> Cost = (vector pvec) => {
        // copy pvec into net.p
        for (int i=0;i<pvec.size;i++) net.p[i] = pvec[i];

        double S = 0.0;
        for (int j = 0; j < M; ++j) {
            double x = xg[j], w = wg[j];
            double F   = net.response(x);
            double F1  = net.dresponse(x);
            double F2  = net.ddresponse(x);

            double xc = x - c;
            double Y   = yc + ypc*xc + xc*xc * F;
            double Y1  = ypc + 2*xc*F + xc*xc * F1;
            double Y2  = 2*F + 4*xc*F1 + xc*xc * F2;

            double R = Phi(Y2, Y1, Y, x);
            S += w * R * R;
        }
        return S; // ≈ ∫ Φ^2
    };

    // Run your central-diff Newton on the scalar Cost
    var resNC = NewtonCentral(Cost, net.p.copy(), 1e-6, 200);
vector pmin = resNC.xmin;
int iters = resNC.iterations;
for (int i = 0; i < pmin.size; i++) net.p[i] = pmin[i];
return (iters, Cost(net.p));

}
static double gprime(double x)
        {
            return Math.Exp(-x * x) * (-5.0 * Math.Sin(5.0 * x - 1.0) - 2.0 * x * Math.Cos(5.0 * x - 1.0));
        }
static double gsecond(double x) {
    // g''(x) = e^{-x^2} [ (4x^2 - 27) cos(5x-1) + 20 x sin(5x-1) ]
    double ex = Math.Exp(-x*x);
    double c  = Math.Cos(5.0*x - 1.0);
    double s  = Math.Sin(5.0*x - 1.0);
    return ex * ( (4.0*x*x - 27.0) * c + 20.0*x * s );
}

        static double g(double x) => Math.Cos(5.0 * x - 1.0) * Math.Exp(-x * x);

        static void Main(string[] args)
        {
            Console.WriteLine("\nTask A:\nIn this task we created an ANN based on our minimization routine from a previous exercise and trained it to approximate the function g(x)=cos(5x-1)*exp(x^2)\n");

            int N = 41;                            // training samples
            double a = -1.0, b = 1.0;
            double[] xs = new double[N];
            double[] ys = new double[N];
            for (int k = 0; k < N; ++k)
            {
                double t = (double)k / (N - 1);
                double x = a + (b - a) * t;
                xs[k] = x;
                ys[k] = g(x);
            }
            int nHidden = 10;
            var net = new ann(nHidden);

            var (iters, mse) = net.train(xs, ys, maxIter: 200, tol: 1e-10);

            Console.WriteLine($"Hidden neurons: {nHidden}");
            Console.WriteLine($"Training samples: {N}");
            Console.WriteLine($"Iterations: {iters}");
            Console.WriteLine($"Final training MSE: {mse:G6}\n");

            Console.WriteLine(" x\t\tg(x)\t\tF_p(x)\t\terror");
            Console.WriteLine("------------------------------------------------------------");
            int T = 11;
            for (int k = 0; k < T; ++k)
            {
                double t = (double)k / (T - 1);
                double x = a + (b - a) * t;
                double gy = g(x);
                double Fy = net.response(x);
                Console.WriteLine($"{x,6:G3}\t{gy,12:G6}\t{Fy,12:G6}\t{(Fy - gy),12:G6}");
            }
            Console.WriteLine();

            using (var w = new System.IO.StreamWriter("ann_data.dat"))
            {
                for (int k = 0; k < N; ++k)
                    w.WriteLine($"{xs[k].ToString("G17", CultureInfo.InvariantCulture)}\t{ys[k].ToString("G17", CultureInfo.InvariantCulture)}");
            }

            using (var w = new System.IO.StreamWriter("ann_fit.dat"))
            {
                int M = 400;
                for (int i = 0; i < M; ++i)
                {
                    double t = (double)i / (M - 1);
                    double x = a + (b - a) * t;
                    double Fy = net.response(x);
                    w.WriteLine($"{x.ToString("G17", CultureInfo.InvariantCulture)}\t{Fy.ToString("G17", CultureInfo.InvariantCulture)}");
                }
            }

            using (var w = new System.IO.StreamWriter("ann_plot.gp"))
            {
                w.WriteLine("set term pngcairo size 900,600");
                w.WriteLine("set out 'ann_fit.png'");
                w.WriteLine("set title 'Task A: ANN fit to g(x)=cos(5x-1)exp(-x^2) with Gaussian wavelets'");
                w.WriteLine("set xlabel 'x'");
                w.WriteLine("set ylabel 'y'");
                w.WriteLine("set key left top");
                w.WriteLine("plot \\");
                w.WriteLine("  'ann_data.dat' using 1:2 with points pt 7 ps 1.2 title 'data', \\");
                w.WriteLine("  'ann_fit.dat'  using 1:2 with lines  lw 2   title 'ANN approximation'");
            }

            Console.WriteLine("To see the full results, check the plots in the files ann_fit.png and ann_residual.png\n");
        

                Console.WriteLine(@"Task B: Modifing the previous exercise such that the network, after training, also can return the first and second derivatives and also the anti-derivative of the approximant to the tabulated function.
    We use F(x) = Σ_i w_i f(z_i),  z_i = (x - a_i)/b_i,  b_i = exp(s_i).
    Since dz_i/dx = 1/b_i and so on and
    For the Gaussian wavelet  f(z) = z·exp(-z^2):
    f'(z)  = (1 - 2 z^2)·exp(-z^2)
    f''(z) = (4 z^3 - 6 z)·exp(-z^2)
    A(z)   = ∫ f(z) dz = -1/2 · exp(-z^2)
    We can use the chain rule to calculate the derivatives and anti-derivatives.

    Therefore (plugging in z_i = (x - a_i)/b_i in the actual method):
    F'(x)  = Σ_i (w_i / b_i) · (1 - 2z_i^2) · exp(-z_i^2)
    F''(x) = Σ_i (w_i / b_i^2) · (4z_i^3 - 6z_i) · exp(-z_i^2)
    ∫_c^x F = -1/2 · Σ_i w_i b_i · [ exp(-z_i(x)^2) - exp(-z_i(c)^2) ]

    Here is a small table where we compare the actual values of the derivatives vs. our approximant's:
    ");

Console.WriteLine(" x\t\tF'(x)\t\tg'(x)\t\terr'\t\tF''(x)\t\tg''(x)\t\terr''");
Console.WriteLine("--------------------------------------------------------------------------------");
for (int k = 0; k < 11; ++k) {
    double t = (double)k / 10.0;
    double x = -1.0 + 2.0 * t;
    double F1 = net.dresponse(x);
    double G1 = gprime(x);
    double F2 = net.ddresponse(x);
    double G2 = gsecond(x);
    Console.WriteLine($"{x,6:G3}\t{F1,12:G6}\t{G1,12:G6}\t{(F1-G1),12:G6}\t{F2,12:G6}\t{G2,12:G6}\t{(F2-G2),12:G6}");
}
Console.WriteLine();
using (var w = new System.IO.StreamWriter("ann_d1.dat")) {
    int M = 400;
    for (int i = 0; i < M; ++i) {
        double t = (double)i / (M - 1);
        double x = -1.0 + 2.0 * t;
        w.WriteLine($"{x.ToString("G17", CultureInfo.InvariantCulture)}\t{net.dresponse(x).ToString("G17", CultureInfo.InvariantCulture)}");
    }
}
using (var w = new System.IO.StreamWriter("ann_d2.dat")) {
    int M = 400;
    for (int i = 0; i < M; ++i) {
        double t = (double)i / (M - 1);
        double x = -1.0 + 2.0 * t;
        w.WriteLine($"{x.ToString("G17", CultureInfo.InvariantCulture)}\t{net.ddresponse(x).ToString("G17", CultureInfo.InvariantCulture)}");
    }
}

Console.WriteLine(" Primitive H(x) = ∫_0^x F(t) dt (ANN only since no simple elementary antiderivative of g):");
for (int k = 0; k < 5; ++k) {
    double[] probe = new double[] { -1.0, -0.5, 0.0, 0.5, 1.0 };
    double x = probe[k];
    Console.WriteLine($"  x={x,5:G3}  H(x)={net.primitive(x, 0.0):G10}");
}
// ---- Task C (ODE) ----
{
    double odeA = 0.0, odeB = Math.PI, odeC = 0.0;
    double yc = 1.0, ypc = 0.0;

    int nHiddenC = 10;
    var netC = new ann(nHiddenC);

    // quick init grid:
    int NinitC = 21;
    double[] xsC = new double[NinitC], ysC = new double[NinitC];
    for (int k = 0; k < NinitC; k++) { xsC[k] = odeA + (odeB - odeA) * k / (NinitC - 1.0); ysC[k] = 0.0; }
    netC.initialize(xsC, ysC);
for (int i = 0; i < nHiddenC; i++)
    netC.p[2 * nHiddenC + i] = 0.1 * Math.Cos((i + 0.5) * Math.PI / nHiddenC);


    var (itersC, JC) = TrainODE(netC, odeA, odeB, M: 201, c: odeC, yc: yc, ypc: ypc);
    Console.WriteLine($"\nTask C:\nHere we use our ANN to approximate a solution to a differential equation. \nWe try it out on the equation y''+y=0 with y(0)=1 and y(0)=0, which actual solution is cos.");
    Console.WriteLine($"Iterations: {itersC}, final integral of residual ≈ {JC:G6}");

    Console.WriteLine(" x\tY(x)\t\tcos(x)\t\terr");
    for (int k = 0; k < 11; k++) {
        double x  = odeA + (odeB - odeA) * k / 10.0;
        double xc = x - odeC;
        double F  = netC.response(x);
        double Y  = yc + ypc * xc + xc * xc * F;
        Console.WriteLine($"{x,6:G3}\t{Y,12:G6}\t{Math.Cos(x),12:G6}\t{(Y - Math.Cos(x)),12:G6}");
    }
}

        }
    }
}
