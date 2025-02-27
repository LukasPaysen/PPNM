using System;
using System.Text;

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
        static void Main(string [] args)
        {

        /*
            if (args.Length > 0 && args[0] == "hello")
    {
        // Run the "hello" part
        
        Console.WriteLine("Hello, this is the first run!");
    }*/
            //Console.WriteLine(args[0]);   
            if (args.Length == 0 || args[0] != "-time")
            {
            // Use a fixed seed for reproducible results.
            Random rnd = new Random(1);

            // ======================================================
            // Test 1: QR Decomposition on a tall matrix A (n > m)
            // ======================================================
            int n = 5, m = 3;
            Console.WriteLine("Testing QR decomposition on a tall matrix A (n > m):");
            matrix A = new matrix(n, m);
            for (int j = 0; j < m; j++)
                for (int i = 0; i < n; i++)
                    A[i, j] = rnd.NextDouble();

            Console.WriteLine("Matrix A:");
            Console.WriteLine(A);

            QR qrTall = new QR(A);

            Console.WriteLine("Matrix Q (from QR decomposition):");
            Console.WriteLine(qrTall.Q);
            matrix QT = new matrix(m,n);
            
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    QT[i,j] = qrTall.Q[j,i];
                }
            }
            Console.WriteLine("Qᵀ:");
            Console.WriteLine(QT);
            Console.WriteLine("Matrix R (from QR decomposition):");
            Console.WriteLine(qrTall.R);

            // Check that R is upper triangular.
            bool isUpper = true;
            for (int i = 1; i < m; i++)
                for (int j = 0; j < i; j++)
                    if (Math.Abs(qrTall.R[i, j]) > 1e-9)
                        isUpper = false;
            Console.WriteLine("R is upper triangular: " + isUpper);


            // Check that Qᵀ*Q = I
            bool qtqIsIdentity = true;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    double dot = 0;
                    for (int k = 0; k < n; k++)
                    {
                        dot += qrTall.Q[k, i] * qrTall.Q[k, j];
                    }
                    if (i == j)
                    {
                        if (Math.Abs(dot - 1) > 1e-9)
                            qtqIsIdentity = false;
                    }
                    else
                    {
                        if (Math.Abs(dot) > 1e-9)
                            qtqIsIdentity = false;
                    }
                }
            }
            Console.WriteLine("Qᵀ * Q:");
            matrix I = new matrix(m,m);
            I = QR.Multiply(QT, qrTall.Q);
            Console.WriteLine(I);
            Console.WriteLine("Qᵀ * Q is identity: " + qtqIsIdentity);
            Console.WriteLine("Q*R:");
            
            matrix testA = new matrix(n,m);
            testA = QR.Multiply(qrTall.Q,qrTall.R);
            Console.WriteLine(testA);


            
            

            // Check that Q*R equals A.
            matrix A_reconstructed = new matrix(n, m);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < m; k++)
                        sum += qrTall.Q[i, k] * qrTall.R[k, j];
                    A_reconstructed[i, j] = sum;
                }
            bool aMatches = true;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    if (Math.Abs(A[i, j] - A_reconstructed[i, j]) > 1e-9)
                        aMatches = false;
            Console.WriteLine("Q * R equals A: " + aMatches);

            // ======================================================
            // Test 2: Solving a linear system using QR decomposition
            // ======================================================
            Console.WriteLine("\nTesting solving linear system with square matrix A:");
            int size = 4;
            matrix A_square = new matrix(size, size);
            vector b = new vector(size);
            for (int j = 0; j < size; j++)
                for (int i = 0; i < size; i++)
                    A_square[i, j] = rnd.NextDouble();
            for (int i = 0; i < size; i++)
                b[i] = rnd.NextDouble();

            Console.WriteLine("Square matrix A:");
            Console.WriteLine(A_square);
            Console.WriteLine("Vector b:");
            Console.WriteLine(b);

            QR qrSquare = new QR(A_square);
            vector x = qrSquare.solve(b);
            Console.WriteLine("Solution vector x (solving A*x = b):");
            Console.WriteLine(x);

            // Check that A * x is (approximately) equal to b.
            vector Ax = new vector(size);
            for (int i = 0; i < size; i++)
            {
                double sum = 0;
                for (int j = 0; j < size; j++)
                    sum += A_square[i, j] * x[j];
                Ax[i] = sum;
            }
            Console.WriteLine("A*x:");
            Console.WriteLine(Ax);
            bool solutionMatches = true;
            for (int i = 0; i < size; i++)
                if (Math.Abs(Ax[i] - b[i]) > 1e-9)
                    solutionMatches = false;
            Console.WriteLine("A*x equals b: " + solutionMatches);


            double detA = qrSquare.det();
            Console.WriteLine("Determinant of A (from R): " + detA);

            Console.WriteLine("\nB.");
            Console.WriteLine("New random square matrix A:");
            matrix A1 = new matrix(m, m);
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    A1[i, j] = rnd.NextDouble();
                }
            }
            Console.WriteLine(A1);
            QR qrTall1 = new QR(A1);
            Console.WriteLine("Matrix Q (from QR decomposition):");
            Console.WriteLine(qrTall1.Q);

            Console.WriteLine("Matrix R (from QR decomposition):");
            Console.WriteLine(qrTall1.R);
            
            Console.WriteLine("The inverse B:");
            matrix iceJJFish =  qrTall1.inverse();
            Console.WriteLine(iceJJFish);
            matrix jjfish = new matrix(m, m);
            jjfish = QR.Multiply(iceJJFish, A1);
            Console.WriteLine("A*B:");
            Console.WriteLine(jjfish);
            Console.WriteLine("So B=A^-1");
            }
            else if (args[0] == "-time")
            {
                Console.Write("Calculating..\t");
                // Default size
                int size = 100;
                // Check for a "-size:" argument to override the default
                foreach (string arg in args)
                {
                    if (arg.StartsWith("-size:"))
                        size = int.Parse(arg.Substring(6));
                }

                // Create a random NxN matrix
                var rnd = new Random(1); // Fixed seed for reproducibility
                matrix A = new matrix(size, size);
                for (int i = 0; i < size; i++)
                {
                    for (int j = 0; j < size; j++)
                    {
                        A[i, j] = rnd.NextDouble();
                    }
                }

                // Time the QR factorization
                //var sw = System.Diagnostics.Stopwatch.StartNew();
                QR qr = new QR(A);
                //sw.Stop();

                // Output the matrix size and the elapsed time in seconds
                //Console.WriteLine($"{size} {sw.Elapsed.TotalSeconds}");
                Console.WriteLine("Done!");
            }



            
            
            
            
            


        }
    }
    
}


