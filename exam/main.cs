// Program.cs
// Compile (example): mcs Program.cs -out:declipping.exe
// Run:               mono declipping.exe
// Output:            declipping.dat

using System;
using System.Text;
using System.IO;
using System.Globalization;
using System.Collections.Generic;

public class vector {
    public double[] data;
    public int size => data.Length;
    public double this[int i] { get => data[i]; set => data[i] = value; }
    public vector(int n){ data = new double[n]; }
    public vector(double[] a){ data = (double[])a.Clone(); }
    public vector copy(){ return new vector((double[])data.Clone()); }
    public override string ToString() => string.Join(", ", data);
}

public class matrix {
    public readonly int size1, size2; // rows, cols
    private double[] data;            // column-major: [i + j*size1]
    public matrix(int n, int m){ size1=n; size2=m; data=new double[n*m]; }
    public double this[int i,int j]{ get=>data[i + j*size1]; set=>data[i + j*size1]=value; }
    public matrix copy(){ var M=new matrix(size1,size2); Array.Copy(data,M.data,data.Length); return M; }
    public override string ToString(){
        var sb=new StringBuilder();
        for(int i=0;i<size1;i++){
            for(int j=0;j<size2;j++) sb.Append(this[i,j].ToString("G6")).Append(j+1<size2? "\t":"");
            sb.AppendLine();
        }
        return sb.ToString();
    }
    // Matrix * vector
    public vector Mul(vector v){
        if(size2!=v.size) throw new Exception("matvec: dim mismatch");
        var r=new vector(size1);
        for(int j=0;j<size2;j++){
            double vj=v[j];
            for(int i=0;i<size1;i++) r[i]+=this[i,j]*vj;
        }
        return r;
    }
    // Extract submatrix of selected columns
    public matrix SubColumns(List<int> cols){
        var A=new matrix(size1, cols.Count);
        for(int k=0;k<cols.Count;k++){
            int j=cols[k];
            for(int i=0;i<size1;i++) A[i,k]=this[i,j];
        }
        return A;
    }
}

public class QR {
    public matrix Q; // n x m, orthonormal columns
    public matrix R; // m x m, upper triangular
    public QR(matrix A){
        int n=A.size1, m=A.size2;
        if(n<m) throw new Exception("QR: need n>=m (tall or square).");
        Q=A.copy();
        R=new matrix(m,m);
        for(int j=0;j<m;j++){
            // v_j currently in Q[:,j]; remove projections on previous q_i
            for(int i=0;i<j;i++){
                double rij=0;
                for(int r=0;r<n;r++) rij += Q[r,i]*Q[r,j];
                R[i,j]=rij;
                for(int r=0;r<n;r++) Q[r,j] -= Q[r,i]*rij;
            }
            // now normalize
            double norm=0;
            for(int r=0;r<n;r++) norm += Q[r,j]*Q[r,j];
            norm = Math.Sqrt(norm);
            if(norm==0) throw new Exception("QR: rank deficiency detected.");
            R[j,j]=norm;
            for(int r=0;r<n;r++) Q[r,j] /= norm;
        }
    }
    public vector solve(vector b){
        // Solve min ||A x - b|| via QR:
        // y = Q^T b (length m), then R x = y
        int n=Q.size1, m=Q.size2;
        var y=new vector(m);
        for(int j=0;j<m;j++){
            double s=0;
            for(int i=0;i<n;i++) s+=Q[i,j]*b[i];
            y[j]=s;
        }
        var x=new vector(m);
        for(int i=m-1;i>=0;i--){
            double s=y[i];
            for(int k=i+1;k<m;k++) s-=R[i,k]*x[k];
            x[i]=s/R[i,i];
        }
        return x;
    }
}

public static class Declipping {

    // Build D (N x N) per lecture notes:
    // rows 0,1: forward 3rd diff [-1, 3, -3, 1] at columns i..i+3
    // rows 2..N-3: centered 3rd derivative [-1/2, 1, 0, -1, 1/2] at columns i-2..i+2
    // rows N-2,N-1: backward 3rd diff [-1, 3, -3, 1] at columns i-3..i
    // Note: scale factors like 1/h^3 are omitted (they don't change argmin).
    public static matrix BuildThirdDerivativeD(int N){
        var D=new matrix(N,N);
        if(N<5) throw new Exception("N must be >= 5 to build 3rd-derivative stencil.");
        // first two rows: forward difference at i=0 and i=1
        for(int i=0;i<2;i++){
            D[i, i+0] += -1.0;
            D[i, i+1] +=  3.0;
            D[i, i+2] += -3.0;
            D[i, i+3] +=  1.0;
        }
        // middle rows: centered difference
        for(int i=2;i<=N-3;i++){
            int j=i-2;
            D[i, j+0] += -0.5;
            D[i, j+1] +=  1.0;
            D[i, j+2] +=  0.0;
            D[i, j+3] += -1.0;
            D[i, j+4] +=  0.5;
        }
        // last two rows: backward difference at i=N-2 and i=N-1
        for(int i=N-2;i<=N-1;i++){
            int j=i-3;
            D[i, j+0] += -1.0;
            D[i, j+1] +=  3.0;
            D[i, j+2] += -3.0;
            D[i, j+3] +=  1.0;
        }
        return D;
    }

    // declip: solve z = argmin || D( y_tilde + M z ) ||, via A z ~= -b, A=D*M, b=D*y_tilde
    // Returns the reconstructed full signal x.
    public static vector Declip(vector y, double yMin, double yMax, bool enforceConsistency=true){
        int N=y.size;
        // 1) find clipped indices in increasing order
        const double eps=1e-12;
        var clipIdx=new List<int>();
        var clipSign=new List<int>(); // +1 for high, -1 for low
        for(int i=0;i<N;i++){
            if(Math.Abs(y[i]-yMax) <= eps){ clipIdx.Add(i); clipSign.Add(+1); }
            else if(Math.Abs(y[i]-yMin) <= eps){ clipIdx.Add(i); clipSign.Add(-1); }
        }
        int n=clipIdx.Count;
        if(n==0) return y.copy(); // nothing to do

        // 2) y_tilde: zeros at clipped positions
        var ytil=y.copy();
        for(int k=0;k<n;k++) ytil[clipIdx[k]] = 0.0;

        // 3) D and b = D * y_tilde
        var D = BuildThirdDerivativeD(N);
        var b = D.Mul(ytil); // size N

        // 4) A = D * M  == columns of D at clipped positions (N x n)
        var A = D.SubColumns(clipIdx);

        // 5) Solve min ||A z + b||, i.e., A z ≈ -b using QR
        //    Do QR on A (N x n), then solve for -b.
        var rhs = b.copy();
        for(int i=0;i<N;i++) rhs[i] = -rhs[i];

        var qr = new QR(A);
        var z  = qr.solve(rhs); // length n

        // 6) Assemble x = y_tilde + M z  (copy y and replace clipped spots)
        var x = y.copy();
        for(int k=0;k<n;k++){
            x[clipIdx[k]] = z[k];
            if(enforceConsistency){
                if(clipSign[k] > 0 && x[clipIdx[k]] < yMax) x[clipIdx[k]] = yMax;
                if(clipSign[k] < 0 && x[clipIdx[k]] > yMin) x[clipIdx[k]] = yMin;
            }
        }
        return x;
    }
}

public class Program {
    public static void Main(){
        // Synthetic test: sine wave clipped at ±0.8; write declipping.dat
        Console.WriteLine("Exam!\nTask 1: \nIn the first part we simply created the function declip that takes a clipped signal, the limits y_min and y_max as inputs and returns the reconstructed signal");

        Console.WriteLine("\nTask 2: ");

        
        
        int N = 1000;
        double yMin = -0.8, yMax = 0.8;
        double cycles = 5.0;   // number of sine periods over [0,1]
        var t = new vector(N);
        var xTrue = new vector(N);
        var yClip = new vector(N);

        for(int i=0;i<N;i++){
            t[i] = (double)i/(N-1);
            xTrue[i] = Math.Sin(2*Math.PI*cycles*t[i]);
            yClip[i] = Math.Max(yMin, Math.Min(yMax, xTrue[i]));
        }

        var xDeclipped = Declipping.Declip(yClip, yMin, yMax, enforceConsistency:true);

        // Write .dat for gnuplot: index  t  x_true  y_clipped  x_declipped  clipped_flag
        string outFile="declipping.dat";
        using(var sw=new StreamWriter(outFile,false,Encoding.UTF8)){
            sw.WriteLine("# i\t t\t x_true\t y_clipped\t x_declipped\t clipped_flag");
            for(int i=0;i<N;i++){
                int clipped = (yClip[i]==yMin || yClip[i]==yMax) ? 1 : 0;
                sw.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                    i, t[i], xTrue[i], yClip[i], xDeclipped[i], clipped));
            }
        }

        // Simple console summary
        int nclip=0; for(int i=0;i<N;i++) if(yClip[i]==yMin || yClip[i]==yMax) nclip++;
        Console.WriteLine("Generated {0} samples; clipped points: {1}", N, nclip);
        Console.WriteLine("\nIn the second part we tested the implementation on the sine wave clipped at -.8 and .8. \nThe result can easily be seen in the file declipping.png\nThe true and declipped signal are indistiguishable in the top plot, which suggests that the declip method performs quite well. \nIn the bottom plot (of the errors) the differences are more apparent, although never very large (|error| always below 0.0002)");

        Console.WriteLine("\nTask 3:");
        // ===== Complex test: chirp + AM + harmonics + transient burst =====
{
    int N2 = 2000;
    double yMin2 = -0.7, yMax2 = 0.7;   // renamed to avoid shadowing

    var t2       = new vector(N2);
    var xTrue2   = new vector(N2);
    var yClip2   = new vector(N2);

    // Ingredients
    double f0 = 3.0;      // starting freq in cycles over [0,1]
    double df = 20.0;     // additional cycles by the end (chirp rate)
    double fH = 25.0;     // harmonic component frequency
    double fL = 1.0;      // slow baseline wobble
    double burstCenter = 0.65;
    double burstSigma  = 0.015; // narrow burst

    for(int i=0;i<N2;i++){
        double tau = 0.6*(double)i/(N2-1);   // renamed from `t` to avoid shadowing
        t2[i] = tau;

        // Chirp: sin(2π [f0 τ + 0.5*df τ^2])
        double chirp = Math.Sin(2*Math.PI*(f0*tau + 0.5*df*tau*tau));

        // Harmonic
        double harmonic = 0.3 * Math.Sin(2*Math.PI*fH*tau);

        // AM envelope (0.2..1.0)
        double env = 0.6 + 0.4 * Math.Cos(2*Math.PI*0.5*tau);

        // Slow baseline wobble
        double baseline = 0.1 * Math.Cos(2*Math.PI*fL*tau);

        // Localized burst (Gaussian window * high-freq sine)
        double gauss = Math.Exp(-0.5*Math.Pow((tau - burstCenter)/burstSigma, 2));
        double burst = 0.3 * gauss * Math.Sin(2*Math.PI*50.0*tau);

        double x = env*(chirp + harmonic) + baseline + burst;
        xTrue2[i] = x;

        // Clip
        yClip2[i] = Math.Max(yMin2, Math.Min(yMax2, x));
    }

    var xDeclipped2 = Declipping.Declip(yClip2, yMin2, yMax2, enforceConsistency:true);

    // Write .dat: i  t  x_true  y_clipped  x_declipped  clipped_flag
    string outFile2="declipping_complex.dat";
    using(var sw=new StreamWriter(outFile2,false,Encoding.UTF8)){
        sw.WriteLine("# i\t t\t x_true\t y_clipped\t x_declipped\t clipped_flag");
        for(int i=0;i<N2;i++){
            int clipped = (yClip2[i]==yMin2 || yClip2[i]==yMax2) ? 1 : 0;
            sw.WriteLine(string.Format(System.Globalization.CultureInfo.InvariantCulture,
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                i, t2[i], xTrue2[i], yClip2[i], xDeclipped2[i], clipped));
        }
    }

    int nclip2=0; for(int i=0;i<N2;i++) if(yClip2[i]==yMin2 || yClip2[i]==yMax2) nclip2++;
    Console.WriteLine("Complex test -> N={0}; clipped points: {1}", N2, nclip2);
    Console.WriteLine("For task 3 we built a tougher signal — a chirp (rising pitch) with a bit of harmonics, a slow amplitude wobble, and a short burst—then clipped at ±0.7 and ran the same declipper, saving/plotting the results in declipping_complex.png");
    Console.WriteLine("We see that the declipper still works well overall—the declipped curve basically sits on top of the true one.The only big miss is ~0.06–0.09 s: that’s a long saturated block with no data, so the method just smooths across it and overshoots the real peak.");
}

    }
}
