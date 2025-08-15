using System;
using System.Collections.Generic;
    using System.IO;


namespace JacobiEVD
{
    class Program
    {
        static void Main(string[] args)
        {
            // Task A: Print header and initial matrix
            Console.WriteLine("Task A:");
            Console.WriteLine("1.");
            Console.WriteLine("Generating random symmetric 4x4 matrix A:");
            Console.WriteLine("A:");

            int size = 4;
            double[,] A = GenerateRandomSymmetricMatrix(size);
            PrintMatrix(A);
            double[,] A1 = CopyMatrix(A);

            // initialize V ← I
            double[,] V = IdentityMatrix(size);

            // run cyclic Jacobi to diagonalize A into A (now D) and accumulate V
            Cyclic(A, V);






            // compute V*D*Vᵀ
            var VDVt = Multiply(Multiply(V, A), Transpose(V));
            var VtAV = Multiply(Multiply(Transpose(V), A1), V);


            // Proving the implementation works by checking that VTAV==D, VDVT==A, VTV==1, VVT==1.
            Console.WriteLine("\n Proving the implementation works by checking that V^TAV==D, VDV^T==A, V^TV==1, VV^T==1. \n We do that by showing that their differences equal the zero matrix. \n We start by simply printing D and V.");
            Console.WriteLine("D:");
            PrintMatrix(A);
            Console.WriteLine("V:");
            PrintMatrix(V);
            Console.WriteLine("\n V^TAV:");
            PrintMatrix(VtAV);

            double[,] C = SubtractMatrices(VtAV, A);
            Console.WriteLine("\nV^TAV - D:");
            PrintMatrix(C);
            Console.WriteLine("\nSo V^TAV==D. We now do the same for the rest - printing each side of the expressions and subtracting.");
            Console.WriteLine("VDV^T:");

            PrintMatrix(VDVt);
            Console.WriteLine("\nA:");

            PrintMatrix(A1);

            Console.WriteLine("\nVDV^T - A:");
            double[,] C1 = SubtractMatrices(VDVt, A1);
            PrintMatrix(C1);

            Console.WriteLine("\nV:");
            PrintMatrix(V);
            Console.WriteLine("\nV^T");
            PrintMatrix(Transpose(V));

            Console.WriteLine("\nV^TV:");
            PrintMatrix(Multiply(Transpose(V), V));

            Console.WriteLine("\nFor completeness we still calculate the difference, V^TV minus the identity matrix. \nV^TV-1:");
            PrintMatrix(SubtractMatrices(Multiply(Transpose(V), V), IdentityMatrix(size)));

            Console.WriteLine("\nVV^T:");
            PrintMatrix(Multiply(V, Transpose(V)));
            Console.WriteLine("\nVV^T - 1:");
            PrintMatrix(SubtractMatrices(Multiply(V, Transpose(V)), IdentityMatrix(size)));



            //Task b
            Console.WriteLine("\nTask B (Numerical Calculation):");
            // read command-line options: e.g. mono main.exe -rmax 10 -dr 0.3
            double rmax = 10.0, dr = 0.3;
            for (int i = 0; i < args.Length; i++)
            {
                if (args[i] == "-rmax") rmax = double.Parse(args[++i]);
                else if (args[i] == "-dr") dr = double.Parse(args[++i]);
            }
            int npoints = (int)(rmax / dr) - 1;

            //  r-grid: r[i]=dr*(i+1), i=0..npoints-1
            double[] r = new double[npoints];
            for (int i = 0; i < npoints; i++)
                r[i] = dr * (i + 1);

            // allocate H and V (V starts as identity)
            double[,] H = new double[npoints, npoints];
            double[,] V1 = IdentityMatrix(npoints);

            // kinetic part:  K = -½·(second-derivative matrix)/dr²
            double diag = -2.0 * (-0.5 / (dr * dr));
            double off = 1.0 * (-0.5 / (dr * dr));
            for (int i = 0; i < npoints - 1; i++)
            {
                H[i, i] = diag;
                H[i, i + 1] = off;
                H[i + 1, i] = off;
            }
            H[npoints - 1, npoints - 1] = diag;

            // add potential W:  W[ii] = -1/r[i]
            for (int i = 0; i < npoints; i++)
            {
                H[i, i] += -1.0 / r[i];
            }

            // diagonalize H in-place: after this H[i,i]=ε_i and V[,i] is the i-th eigenvector
            Cyclic(H, V1);

            // collect eigenpairs into a list so we can sort by energy
            var eigenpairs = new List<(double energy, int idx)>();
            for (int i = 0; i < npoints; i++)
            {
                eigenpairs.Add((H[i, i], i));
            }
            // sort by increasing energy (most negative first)
            eigenpairs.Sort((a, b) => a.energy.CompareTo(b.energy));

            // print out the first few eigenvalues and eigenvectors
            int numToShow = Math.Min(5, npoints);
            Console.WriteLine("\nLowest {0} eigenvalues and corresponding eigenvectors (from rmax = {1} and dr = {2}):", numToShow, rmax, dr);
            for (int k = 0; k < numToShow; k++)
            {
                var (energy, idx) = eigenpairs[k];
                Console.WriteLine($"\n  k={k}  ε_{k} = {energy:F6}  Ha   (grid-index = {idx})");
                Console.Write("    f^{(k)}_i = 1/sqrt(dr) * V[i,idx]: ");
                for (int i = 0; i < Math.Min(10, npoints); i++)
                {      // show only first 10 components
                    double fval = V1[i, idx] / Math.Sqrt(dr);
                    Console.Write($"{fval:F4} ");
                }
                Console.WriteLine("…");
            }

            Console.WriteLine("\n(Convergence): \nWe investigate convergence by plotting the lowest eigenvalue obtained against different values of rmax and dr. See the files \"epsilon0_vs_dr.png\" and \"epsilon0_vs_rmax.png\"\n the analytical result for ε_0 is -0.5.");


            // --- begin convergence‐in‐dr block ---

            // fix rmax, write 100 (dr, ε0) pairs to a tab-separated file
            double fixedRmax = 10.0;
            int points = 20;
            using (var w = new StreamWriter("conv_dr_vs_eps0.dat"))
            {
                for (int j = 0; j < points; j++)
                {
                    // dr runs from 0.01 to 1.00
                    double drj = 0.01 + j * (1.0 - 0.01) / (points - 1);
                    int npts = (int)(fixedRmax / drj) - 1;

                    // build H for this dr
                    double[] rj = new double[npts];
                    for (int i = 0; i < npts; i++) rj[i] = drj * (i + 1);

                    var Hj = new double[npts, npts];
                    double diagj = -2.0 * (-0.5 / (drj * drj));
                    double offj = 1.0 * (-0.5 / (drj * drj));
                    for (int i = 0; i < npts - 1; i++)
                    {
                        Hj[i, i] = diagj;
                        Hj[i, i + 1] = offj;
                        Hj[i + 1, i] = offj;
                    }
                    Hj[npts - 1, npts - 1] = diagj;
                    for (int i = 0; i < npts; i++)
                        Hj[i, i] += -1.0 / rj[i];

                    // diagonalize
                    double[,] Vtmp = IdentityMatrix(npts);
                    Cyclic(Hj, Vtmp);

                    // find lowest eigenvalue (smallest H[i,i])
                    double eps0 = double.PositiveInfinity;
                    for (int i = 0; i < npts; i++)
                        if (Hj[i, i] < eps0) eps0 = Hj[i, i];

                    // write dr and ε0
                    w.WriteLine($"{drj:F6}\t{eps0:F6}");
                }
            }
            // --- end convergence‐in‐dr block ---
            // --- convergence‐in‐rmax block ---
            using (var w2 = new StreamWriter("conv_rmax_vs_eps0.dat"))
            {
                double fixedDr = 0.1;        // choose a reasonable grid‐spacing
                int samples = 20;            // number of rmax values

                for (int j = 0; j < samples; j++)
                {
                    // sweep rmax from 5.0 to 20.0
                    double rmaxj = 3.0 + j * (20.0 - 3.0) / (samples - 1);
                    int npts = (int)(rmaxj / fixedDr) - 1;

                    // build grid
                    double[] rj = new double[npts];
                    for (int i = 0; i < npts; i++) rj[i] = fixedDr * (i + 1);

                    // assemble H
                    double[,] Hj = new double[npts, npts];
                    double diaggg = -2.0 * (-0.5 / (fixedDr * fixedDr));
                    double offff = 1.0 * (-0.5 / (fixedDr * fixedDr));
                    for (int i = 0; i < npts - 1; i++)
                    {
                        Hj[i, i] = diaggg;
                        Hj[i, i + 1] = offff;
                        Hj[i + 1, i] = offff;
                    }
                    Hj[npts - 1, npts - 1] = diaggg;
                    for (int i = 0; i < npts; i++)
                        Hj[i, i] += -1.0 / rj[i];

                    // diagonalize
                    var Vtmp = IdentityMatrix(npts);
                    Cyclic(Hj, Vtmp);

                    // extract the lowest eigenvalue
                    double eps0 = double.PositiveInfinity;
                    for (int i = 0; i < npts; i++)
                        if (Hj[i, i] < eps0) eps0 = Hj[i, i];

                    // write rmax and ε0
                    w2.WriteLine($"{rmaxj:F6}\t{eps0:F6}");
                }
            }
            // --- end block ---

            // --- write 3 numerics + 1 analytic to wave_functions.dat ---
            using (var w = new StreamWriter("wave_functions.dat"))
            {
                w.WriteLine("# r\tf_num0\tf_num1\tf_num2\tf_exact0");
                for (int i = 0; i < npoints; i++)
                {
                    double ri = dr * (i + 1);
                    // find grid‐indices of the three lowest states
                    int idx0 = eigenpairs[0].idx;
                    int idx1 = eigenpairs[1].idx;
                    int idx2 = eigenpairs[2].idx;

                    // numerical reduced wave-functions
                    double f0_num = V1[i, idx0] / Math.Sqrt(dr);
                    double f1_num = V1[i, idx1] / Math.Sqrt(dr);
                    double f2_num = V1[i, idx2] / Math.Sqrt(dr);

                    // analytic ground-state: f0(r)=2 r e^{-r}
                    double f0_exact = 2.0 * ri * Math.Exp(-ri);
                    double f1_exact = (1.0 / Math.Sqrt(2.0)) * (2 * ri - ri * ri) * Math.Exp(-ri / 2.0);

                    w.WriteLine($"{ri:F6}\t{f0_num:F6}\t{f1_num:F6}\t{f2_num:F6}\t{f0_exact:F6}\t{f1_exact:F6}");
                }
            }
            Console.WriteLine("(Wave-functions):");
            Console.WriteLine("The lowest eigenfunctions are plotted and compared to analytical results in \"wavefunctions.png\"");
            













        }
        static double[,] SubtractMatrices(double[,] A, double[,] B)
        {
            int rows = A.GetLength(0), cols = A.GetLength(1);
            double[,] result = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    result[i, j] = A[i, j] - B[i, j];
            return result;
        }

        static double[,] CopyMatrix(double[,] A)
        {
            int rows = A.GetLength(0), cols = A.GetLength(1);
            double[,] copy = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    copy[i, j] = A[i, j];
            return copy;
        }

        static double[,] GenerateRandomSymmetricMatrix(int size)
        {
            Random rand = new Random();
            double[,] matrix = new double[size, size];

            // Fill only the upper triangle and mirror it to the lower triangle
            for (int i = 0; i < size; i++)
            {
                for (int j = i; j < size; j++)
                {
                    double value = rand.NextDouble() * 10 - 5; // Random number in range [-5, 5)
                    matrix[i, j] = value;
                    matrix[j, i] = value;
                }
            }

            return matrix;
        }
        static void PrintMatrix(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    Console.Write($"{matrix[i, j]:F2}\t");
                }
                Console.WriteLine();
            }
        }
        // multiplies A ← A·J(p,q,θ)
        static void timesJ(double[,] A, int p, int q, double theta)
        {
            double c = Math.Cos(theta), s = Math.Sin(theta);
            int n = A.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                double aip = A[i, p], aiq = A[i, q];
                A[i, p] = c * aip - s * aiq;
                A[i, q] = s * aip + c * aiq;
            }
        }

        // multiplies A ← Jᵀ(p,q,θ)·A
        static void Jtimes(double[,] A, int p, int q, double theta)
        {
            double c = Math.Cos(theta), s = Math.Sin(theta);
            int n = A.GetLength(0);
            for (int j = 0; j < n; j++)
            {
                double apj = A[p, j], aqj = A[q, j];
                A[p, j] = c * apj + s * aqj;
                A[q, j] = -s * apj + c * aqj;
            }
        }

        // Build an n×n identity matrix
        static double[,] IdentityMatrix(int n)
        {
            var I = new double[n, n];
            for (int i = 0; i < n; i++) I[i, i] = 1.0;
            return I;
        }

        // The cyclic Jacobi sweep: diagonalizes A into D (in-place) and accumulates rotations in V
        static void Cyclic(double[,] A, double[,] V)
        {
            int n = A.GetLength(0);
            bool changed;
            do
            {
                changed = false;
                for (int p = 0; p < n - 1; p++)
                {
                    for (int q = p + 1; q < n; q++)
                    {
                        double apq = A[p, q], app = A[p, p], aqq = A[q, q];
                        double theta = 0.5 * Math.Atan2(2 * apq, aqq - app);
                        double c = Math.Cos(theta), s = Math.Sin(theta);
                        // compute what the new diagonals would be
                        double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
                        double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;
                        if (new_app != app || new_aqq != aqq)
                        {
                            changed = true;
                            timesJ(A, p, q, theta);      // A ← A·J
                            Jtimes(A, p, q, -theta);      // A ← Jᵀ·A
                            timesJ(V, p, q, theta);      // V ← V·J
                        }
                    }
                }
            } while (changed);
        }


        /// <summary>
        /// In-place multiplies the given matrix A by the Jacobi rotation J(p,q,theta) from the right: A <- A * J
        /// </summary>
        /// <param name="A">Square matrix A (n x n)</param>
        /// <param name="p">First rotation index</param>
        /// <param name="q">Second rotation index</param>
        /// <param name="theta">Rotation angle</param>


        // TODO: Implement Jtimes for left-multiplication J * A
        // TODO: Add cyclic sweep function to perform full Jacobi eigenvalue algorithm
        // helper: multiply two n×n matrices
        static double[,] Multiply(double[,] A, double[,] B)
        {
            int n = A.GetLength(0);
            var C = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    for (int k = 0; k < n; k++)
                        C[i, j] += A[i, k] * B[k, j];
            return C;
        }

        // helper: transpose an n×n matrix
        static double[,] Transpose(double[,] A)
        {
            int n = A.GetLength(0);
            var T = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    T[j, i] = A[i, j];
            return T;
        }

    }
    
}
