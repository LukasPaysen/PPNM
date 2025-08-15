    using System;
    using System.Text;
    using System.IO;
    using System.Globalization;

    namespace QRDecompositionExample
    {
        // A simple vector class
        public class vector
        {
            public double[] data;
            public int size => data.Length;

            public double this[int i]
            {
                get => data[i];
                set => data[i] = value;
            }

            public vector(int n)
            {
                data = new double[n];
            }

            public override string ToString()
            {
                return string.Join(", ", data);
            }
        }

        // A simple matrix class (stored in column-major order)
        public class matrix
        {
            public readonly int size1, size2; // number of rows and columns
            private double[] data;  // internal storage

            public matrix(int n, int m)
            {
                size1 = n;
                size2 = m;
                data = new double[n * m];
            }

            // Indexer: element at row i, column j
            public double this[int i, int j]
            {
                get => data[i + j * size1];
                set => data[i + j * size1] = value;
            }

            // Create a copy of the matrix
            public matrix copy()
            {
                matrix M = new matrix(size1, size2);
                Array.Copy(data, M.data, data.Length);
                return M;
            }

            public override string ToString()
            {
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < size1; i++)
                {
                    for (int j = 0; j < size2; j++)
                    {
                        sb.Append(this[i, j].ToString("F3") + "\t");
                    }
                    sb.AppendLine();
                }
                return sb.ToString();
            }
        }

        // The QR class implements modified Gram-Schmidt decomposition,
        // solving linear systems and computing determinants.
        public class QR
        {
            public matrix Q; // n x m matrix with orthonormal columns
            public matrix R; // m x m upper-triangular matrix

            // Constructor: perform the QR-decomposition on matrix A (n x m with n>=m)
            public QR(matrix A)
            {
                int n = A.size1;
                int m = A.size2;
                // Copy A into Q (which will be modified)
                Q = A.copy();
                R = new matrix(m, m); // R will be m x m

                // Modified Gram-Schmidt
                for (int j = 0; j < m; j++)
                {
                    // Compute the norm of column j of Q
                    double norm = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        norm += Q[i, j] * Q[i, j];
                    }
                    norm = Math.Sqrt(norm);
                    R[j, j] = norm;

                    // Normalize the j-th column of Q
                    for (int i = 0; i < n; i++)
                    {
                        Q[i, j] /= norm;
                    }

                    // Orthogonalize the remaining columns against the j-th column
                    for (int k = j + 1; k < m; k++)
                    {
                        double dot = 0.0;
                        for (int i = 0; i < n; i++)
                        {
                            dot += Q[i, j] * Q[i, k];
                        }
                        R[j, k] = dot;
                        for (int i = 0; i < n; i++)
                        {
                            Q[i, k] -= Q[i, j] * dot;
                        }
                    }
                }
            }
            public matrix inverse()
            {
                int n = Q.size1; // Assuming A is square, so Q is n x n
                matrix inv = new matrix(n, n);

                // Compute each column of the inverse by solving A * x = e,
                // where e is a column of the identity matrix.
                for (int i = 0; i < n; i++)
                {
                    vector e = new vector(n);
                    e[i] = 1.0; // Set the i-th entry to 1 (identity matrix column)
                    vector x = solve(e); // Use th2e already defined solve method
                    // Place the solution vector x into the i-th column of the inverse matrix
                    for (int j = 0; j < n; j++)
                    {
                        inv[j, i] = x[j];
                    }
                }
                return inv;
            }

            // Solve the system Ax = b, where A = Q*R. Since Q is orthogonal,
            // we solve R*x = Qᵀ*b by back-substitution.
            public static matrix Multiply(matrix A, matrix B)
            {
                if (A.size2 != B.size1)
                    throw new Exception("Incompatible dimensions for multiplication.");

                matrix C = new matrix(A.size1, B.size2);
                for (int i = 0; i < A.size1; i++)
                    for (int j = 0; j < B.size2; j++)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < A.size2; k++)
                            sum += A[i, k] * B[k, j];
                        C[i, j] = sum;
                    }
                return C;
            }
            public vector solve(vector b)
            {
                int n = Q.size1;
                int m = Q.size2; // A is n x m (here we assume A is square when solving)
                vector y = new vector(m);

                // Compute y = Qᵀ * b
                for (int j = 0; j < m; j++)
                {
                    double sum = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        sum += Q[i, j] * b[i];
                    }
                    y[j] = sum;
                }

                // Back-substitution to solve R * x = y
                vector x = new vector(m);
                for (int j = m - 1; j >= 0; j--)
                {
                    double sum = y[j];
                    for (int k = j + 1; k < m; k++)
                    {
                        sum -= R[j, k] * x[k];
                    }
                    x[j] = sum / R[j, j];
                }
                return x;
            }

            // Determinant of R (an upper-triangular matrix) is the product of its diagonal elements.
            // (For square A, this is equal to det(A).)
            public double det()
            {
                int m = R.size1;
                double product = 1.0;
                for (int i = 0; i < m; i++)
                {
                    product *= R[i, i];
                }
                return product;
            }
        }
        

        class Program
        {

        // Put these methods in your Program class (or a new static class).

        // Least-squares using weighted design matrix built from fs.
        // Expects:
        //   fs: array of basis functions f_k(x)
        //   x: vector of x_i (times)
        //   y: vector of observations (here: ln(y_i))
        //   dy: vector of uncertainties of y (here: sigma_ln = delta_y / y)
        // Returns (c, Cov) where c is vector of fitted coefficients and Cov is covariance matrix.
        static (vector, matrix) lsfit(Func<double,double>[] fs, vector x, vector y, vector dy)
        {
            int n = x.size;
            int m = fs.Length;
            var A = new matrix(n, m);
            var b = new vector(n);

            // Build weighted system: divide each row by dy[i]
            for (int i = 0; i < n; i++)
            {
                double sigma = dy[i];
                b[i] = y[i] / sigma;
                for (int k = 0; k < m; k++)
                    A[i, k] = fs[k](x[i]) / sigma;
            }

            // QR and solve
            QR qr = new QR(A);
            vector c = qr.solve(b); // solves R c = Q^T b

            // Extract R (m x m)
            matrix R = qr.R;

            // Invert upper-triangular R by solving R * X = I (column by column)
            matrix Rinv = new matrix(m, m);
            for (int col = 0; col < m; col++)
            {
                // Solve R * xcol = e_col
                var xcol = new vector(m);
                for (int i = m - 1; i >= 0; i--)
                {
                    double s = (i == col) ? 1.0 : 0.0;
                    for (int j = i + 1; j < m; j++)
                        s -= R[i, j] * xcol[j];
                    xcol[i] = s / R[i, i];
                }
                for (int i = 0; i < m; i++) Rinv[i, col] = xcol[i];
            }

            // Cov = Rinv * Rinv^T
            matrix RinvT = Transpose(Rinv);
            matrix Cov = QR.Multiply(Rinv, RinvT);

            return (c, Cov);
        }

        // Helper: transpose a matrix
        static matrix Transpose(matrix M)
        {
            matrix T = new matrix(M.size2, M.size1);
            for (int i = 0; i < M.size1; i++)
                for (int j = 0; j < M.size2; j++)
                    T[j, i] = M[i, j];
            return T;
        }

            // Example runner for the Rutherford ThX data (computes c0, c1, a, lambda, T1/2)
            // Call RunThXFit() from Main() to execute.
            static void RunThXFit()
            {
                // Raw data from the exercise
                double[] t_arr = new double[] { 1, 2, 3, 4, 6, 9, 10, 13, 15 };
                double[] y_arr = new double[] { 117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1 };
                double[] dy_arr = new double[] { 6, 5, 4, 4, 4, 3, 3, 2, 2 }; // approx uncertainties in y

                int n = t_arr.Length;
                var x = new vector(n);
                var ylog = new vector(n);   // ln(y)
                var dylog = new vector(n);  // sigma(ln y) = dy / y

                for (int i = 0; i < n; i++)
                {
                    x[i] = t_arr[i];
                    ylog[i] = Math.Log(y_arr[i]);
                    dylog[i] = dy_arr[i] / y_arr[i];
                }

                // Basis functions for ln(y) = c0 + c1 * t  (we choose f1 = -t so c1 = -lambda)
                Func<double, double>[] fs = new Func<double, double>[]
                {
                t => 1.0,
                t => -t
                };

                // Do the fit
                (vector c, matrix Cov) = lsfit(fs, x, ylog, dylog);

                // Extract results
                double c0 = c[0];
                double c1 = c[1];
                double lnA = c0;
                double lambda = c1;               // because we used -t in basis
                double A = Math.Exp(lnA);
                double halfLife = Math.Log(2.0) / lambda;

                // Errors (std dev) from covariance (diagonal)
                double sigma_c0 = Math.Sqrt(Math.Abs(Cov[0, 0]));
                double sigma_c1 = Math.Sqrt(Math.Abs(Cov[1, 1]));
                double sigma_lnA = sigma_c0;
                double sigma_A = A * sigma_lnA;   // propagate: sigma_a = a * sigma_lnA
                double sigma_lambda = sigma_c1;   // lambda = -c1, same uncertainty
                double sigma_halfLife = (Math.Log(2.0) / (lambda * lambda)) * sigma_lambda; // propagate

                // Print
                Console.WriteLine("Task A:");
                Console.WriteLine("Fit results (ln(y) = c0 - c1 * t):");
                Console.WriteLine($"c0 = {c0:F6}");
                Console.WriteLine($"c1 = {c1:F6}");
                Console.WriteLine();
                Console.WriteLine("Derived (original) parameters:");
                Console.WriteLine($"a = exp(c0) = {A:F6}");
                Console.WriteLine($"lambda = {lambda:F6} ");
                Console.WriteLine($"T1/2 = ln(2)/lambda = {halfLife:F6} days");
                Console.WriteLine("The modern value for 224-Ra is 3.63 days while we find it to be 4.06 days. This turns out to be a difference of 11.9 %");
                Console.WriteLine("To see the plot with the experimental data, error bars, and the best fit, open the file thx_ln_plot.png.\n");
                
                Console.WriteLine("Task b:");
                Console.WriteLine("We modify the lsfit function to also calculate the covariant matrix and the uncertainties of the fitting coefficients. The uncertainties are simply the square roots of the diagonal elements of the covariant matrix. We get:");
                Console.WriteLine($"c0 = {c0:F4} ± {sigma_c0:F4}");
                Console.WriteLine($"c1 = {c1:F4} ± {sigma_c1:F4}");

                Console.WriteLine("We want to also estimate the uncertainty of the half-life of ThX. We do that by error-propagation as the error of the half-life T_½ comes from the error in lambda, sigma_lambda, by the formula: \n sigma_T_½ = |d T_½ / d lambda|*sigma_lambda = |-ln(2)/lambda^2|*sigma_lambda");
                Console.WriteLine($"This gives: \n sigma_T_½ = {sigma_halfLife:F4} days \nBut since our half-life minus the uncertainty is still larger than the modern value: \n {halfLife:F4} - {sigma_halfLife:F4} days = {halfLife - sigma_halfLife} > 3.63 days \nit does not agree with the modern value");

                // Task C: write perturbed-fit curves (minimal snippet)
Console.WriteLine("\nTask C: \nWe plot the the best fit function along with four more graphs where we change the coeffiecients by the uncertainties we calculated in task B. \nSee the result in the file thx_ln_unc_plot.png (+/+ means both coefficients are modified by adding their uncertainties to them and so on)");

// compute 1-sigma uncertainties for coefficients
int m = c.size;
var sigma_c = new vector(m);
for (int i = 0; i < m; i++) sigma_c[i] = Math.Sqrt(Math.Abs(Cov[i, i]));

// grid for smooth curves (same range as fit file)
int npts = 200;
double tminC = 0.0;
double tmaxC = 16.0;
double dt = (tmaxC - tminC) / Math.Max(1, npts - 1);

using (var sw = new StreamWriter("thx_ln_perturbed.dat", false, System.Text.Encoding.UTF8))
{
    // header
    sw.WriteLine("# t\tln_central\tln_pp\tln_pm\tln_mp\tln_mm");
    for (int k = 0; k < npts; k++)
    {
        double tval = tminC + k * dt;

        // central
        double ln_central = c[0] - c[1] * tval; // model: ln y = c0 - c1*t (matches your WriteLnFitFile)

        // four sign combos: (+,+), (+,-), (-,+), (-,-)
        double c0_pp = c[0] + sigma_c[0];
        double c1_pp = c[1] + sigma_c[1];
        double ln_pp = c0_pp - c1_pp * tval;

        double c0_pm = c[0] + sigma_c[0];
        double c1_pm = c[1] - sigma_c[1];
        double ln_pm = c0_pm - c1_pm * tval;

        double c0_mp = c[0] - sigma_c[0];
        double c1_mp = c[1] + sigma_c[1];
        double ln_mp = c0_mp - c1_mp * tval;

        double c0_mm = c[0] - sigma_c[0];
        double c1_mm = c[1] - sigma_c[1];
        double ln_mm = c0_mm - c1_mm * tval;

        sw.WriteLine(
            tval.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
            ln_central.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
            ln_pp.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
            ln_pm.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
            ln_mp.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
            ln_mm.ToString("G17", CultureInfo.InvariantCulture)
        );
    }
}



                // where you have the arrays and fitted vector c:
                string datPath = "thx_ln_data.dat";
    string fitPath = "thx_ln_fit.dat";

    // t_arr, y_arr, dy_arr defined earlier in RunThXFit()
    WriteLnDataFile(datPath, t_arr, y_arr, dy_arr);

    // choose t-range and resolution for the smooth fitted line
    double tmin = 0.0;
    double tmax = 16.0; // or Math.Max(t_arr) + 1
    WriteLnFitFile(fitPath, tmin, tmax, 200, c);



            }

    // Writes the experimental data file: t \t ln(y) \t delta_ln(y)
    static void WriteLnDataFile(string path, double[] t_arr, double[] y_arr, double[] dy_arr)
    {
        using (var sw = new StreamWriter(path, false, System.Text.Encoding.UTF8))
        {
            // header (gnuplot ignores lines starting with #)
            sw.WriteLine("# t\tln(y)\tdelta_ln(y)");
            for (int i = 0; i < t_arr.Length; i++)
            {
                double t = t_arr[i];
                double y = y_arr[i];
                double dy = dy_arr[i];
                double lnY = Math.Log(y);
                double dlnY = dy / y; // delta ln(y) = delta y / y
                sw.WriteLine(
                    t.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
                    lnY.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
                    dlnY.ToString("G17", CultureInfo.InvariantCulture)
                );
            }
        }
    }

    // Writes a dense fit-line file for plotting the fitted ln(y) curve
    // c is the fitted coefficient vector c[0]=c0, c[1]=c1 (for ln y = c0 + c1*t)
    static void WriteLnFitFile(string path, double tmin, double tmax, int npts, vector c)
    {
        using (var sw = new StreamWriter(path, false, System.Text.Encoding.UTF8))
        {
            sw.WriteLine("# t\tln(y_fit)");
            double dt = (tmax - tmin) / Math.Max(1, npts - 1);
            for (int k = 0; k < npts; k++)
            {
                double t = tmin + k * dt;
                double lnfit = c[0] - c[1] * t; // if you used basis [1, -t] then use c[0] + c[1]*(-t) instead
                sw.WriteLine(
                    t.ToString("G17", CultureInfo.InvariantCulture) + "\t" +
                    lnfit.ToString("G17", CultureInfo.InvariantCulture)
                );
            }
        }
    }


            static void Main(string[] args)
            {
                RunThXFit();

            }
        }
        
    }


