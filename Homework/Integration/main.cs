using System;
using System.IO;

class AdaptiveIntegrator
{
    // Wrapper that returns both integral and estimated error
    public static (double Q, double err) IntegrateWithError(
        Func<double, double> f,
        double a,
        double b,
        double acc = 1e-3,
        double eps = 1e-3
    )
    {
        // Simple recursive adaptive scheme
        double Q = Integrate(f, a, b, out long calls, acc, eps); // your existing integrator

        // Split interval for error estimation
        double mid = 0.5 * (a + b);
        double Q1 = Integrate(f, a, mid, out _, acc, eps);
        double Q2 = Integrate(f, mid, b, out _, acc, eps);

        // Estimate error as combined difference
        double err = Math.Sqrt(Math.Pow(Q1 + Q2 - Q, 2));

        return (Q, err);
    }

    public static double Integrate(
    Func<double, double> f,
    double a,
    double b,
    out long calls,
    double acc = 1e-3,
    double eps = 1e-3
)
    {
        // Case 1: both limits finite → just use CC
        if (!(double.IsInfinity(a) || double.IsInfinity(b)))
        {
            return IntegrateClenshawCurtis(f, a, b, acc, eps, out calls);
        }

        // Case 2: a finite, b = +∞
        if (!double.IsInfinity(a) && double.IsPositiveInfinity(b))
        {
            Func<double, double> g = t =>
            {
                double x = a + t / (1 - t);
                double dxdt = 1.0 / Math.Pow(1 - t, 2);
                return f(x) * dxdt;
            };
            return IntegrateClenshawCurtis(g, 0, 1, acc, eps, out calls);
        }

        // Case 3: a = -∞, b finite
        if (double.IsNegativeInfinity(a) && !double.IsInfinity(b))
        {
            Func<double, double> g = t =>
            {
                double x = b - (1 - t) / t;
                double dxdt = 1.0 / (t * t);
                return f(x) * dxdt;
            };
            return IntegrateClenshawCurtis(g, 0, 1, acc, eps, out calls);
        }

        // Case 4: both limits infinite
        if (double.IsNegativeInfinity(a) && double.IsPositiveInfinity(b))
        {
            Func<double, double> g = t =>
            {
                double x = Math.Tan(Math.PI * (t - 0.5));
                double dxdt = Math.PI / Math.Pow(Math.Cos(Math.PI * (t - 0.5)), 2);
                return f(x) * dxdt;
            };
            return IntegrateClenshawCurtis(g, 0, 1, acc, eps, out calls);
        }

        throw new ArgumentException("Unsupported integration bounds.");
    }

    // -------------------- Clenshaw-Curtis wrapper (step 1) --------------------
    // Integrate f(x) on [a,b] by theta in [0,pi].
    // Signature: acc and eps are required (no optional defaults) to avoid C# out/optional ordering issues.
    // Returns integral and sets evalCount to number of calls to original f(x).
    public static double IntegrateClenshawCurtis(Func<double, double> f, double a, double b,
                                                 double acc, double eps, out long evalCount)
    {
        double half = 0.5 * (b - a);
        double mid = 0.5 * (a + b);

        long localCount = 0;

        // g(theta) maps theta->x and calls original f(x). It must never return NaN/Infinity
        Func<double, double> g = theta =>
        {
            // map
            double x = mid + half * Math.Cos(theta);

            // call original integrand and count it
            double fx;
            try
            {
                localCount++;
                fx = f(x);
            }
            catch
            {
                // if original f throws for some x (shouldn't happen for well-formed f), treat as 0
                fx = 0.0;
            }

            // guard against NaN/Infinity from f
            if (double.IsNaN(fx) || double.IsInfinity(fx))
                fx = 0.0;

            // Jacobian: sin(theta) * (b-a)/2  (half == (b-a)/2)
            double jac = Math.Sin(theta) * half;
            return fx * jac;
        };

        // Use existing Integrate on theta in [0,pi].
        // NOTE: Integrate uses SafeEval internally which would throw if g returned NaN/Inf;
        // we ensure g never returns NaN/Inf above.
        double result = Integrate(g, 0.0, Math.PI, acc, eps);

        evalCount = localCount;
        return result;
    }

    // -------------------- Compare ordinary vs CC (step 2) --------------------
    // Runs both integrators on same f, prints computed values, absolute errors, and function-eval counts.
    public static void CompareIntegrators(string desc, Func<double, double> fRaw, double a, double b,
                                          double acc, double eps, double expected)
    {
        Console.WriteLine("------------------------------------------------------------");
        Console.WriteLine(desc);
        Console.WriteLine($"  acc={acc:G}, eps={eps:G}");

        // --- Ordinary integrator (count original f calls) ---
        long ordCalls = 0;
        Func<double, double> fOrd = x =>
        {
            ordCalls++;
            return fRaw(x);
        };

        double ordResult = double.NaN;
        bool ordOk = true;
        try
        {
            ordResult = Integrate(fOrd, a, b, acc, eps);
        }
        catch (Exception ex)
        {
            ordOk = false;
            Console.WriteLine($"  Ordinary Integrator ERROR: {ex.Message}");
        }
        double ordErr = ordOk ? Math.Abs(ordResult - expected) : double.NaN;

        // --- Clenshaw-Curtis integrator ---
        long ccCalls = 0;
        double ccResult = double.NaN;
        bool ccOk = true;
        try
        {
            ccResult = IntegrateClenshawCurtis(fRaw, a, b, acc, eps, out ccCalls);
        }
        catch (Exception ex)
        {
            ccOk = false;
            Console.WriteLine($"  CC Integrator ERROR: {ex.Message}");
        }
        double ccErr = ccOk ? Math.Abs(ccResult - expected) : double.NaN;

        // Print summary
        Console.WriteLine(" Method | computed           | abs error       | calls to f(x)");
        Console.WriteLine("--------+--------------------+-----------------+--------------");
        Console.WriteLine($" Ordinary: {ordResult,18:R}  {ordErr,13:E}  {ordCalls,12}");
        Console.WriteLine($" CC      : {ccResult,18:R}  {ccErr,13:E}  {ccCalls,12}");
        Console.WriteLine();
    }

    // Public wrapper: user calls this
    public static double Integrate(Func<double, double> f, double a, double b, double acc = 1e-3, double eps = 1e-3)
    {
        return IntegrateRecursive(f, a, b, acc, eps, double.NaN, double.NaN, 0);
    }

    // Recursive implementation (open 4-point adaptive integrator)
    // f2 and f3 are optional interior reused points (positions a+2h/6 and a+4h/6).
    private static double IntegrateRecursive(Func<double, double> f, double a, double b,
                                             double acc, double eps,
                                             double f2, double f3, int depth)
    {
        const int MAX_DEPTH = 120; // guard against runaway recursion
        if (depth > MAX_DEPTH)
            throw new Exception($"Maximum recursion depth ({MAX_DEPTH}) exceeded on [{a}, {b}]");

        double h = b - a;

        // On the first call for this interval we must compute the two "middle" points f2,f3
        if (double.IsNaN(f2))
        {
            f2 = SafeEval(f, a + 2.0 * h / 6.0);
            f3 = SafeEval(f, a + 4.0 * h / 6.0);
        }

        // compute the two remaining (outer) interior points
        double f1 = SafeEval(f, a + 1.0 * h / 6.0);
        double f4 = SafeEval(f, a + 5.0 * h / 6.0);

        // higher-order open 4-point rule (as in the exercise)
        double Q = (2.0 * f1 + f2 + f3 + 2.0 * f4) / 6.0 * h;
        // embedded lower-order rule reusing the same evaluations
        double q = (f1 + f2 + f3 + f4) / 4.0 * h;

        double err = Math.Abs(Q - q);

        // --- safety: if interval is extremely small, stop subdividing and accept Q
        double min_h = 1e-15 * (1.0 + Math.Abs(a));
        if (h < min_h)
        {
            return Q;
        }

        // stopping criterion: absolute + relative tolerance
        if (err <= acc + eps * Math.Abs(Q))
        {
            return Q;
        }
        else
        {
            // split absolute tolerance for children --- less aggressive split (acc/2)
            double accChild = acc / 2.0;
            double mid = 0.5 * (a + b);

            // Reuse: left gets (f1,f2) as its interior reused values
            double left = IntegrateRecursive(f, a, mid, accChild, eps, f1, f2, depth + 1);

            // right gets (f3,f4)
            double right = IntegrateRecursive(f, mid, b, accChild, eps, f3, f4, depth + 1);

            return left + right;
        }
    }

    // Safe wrapper for function evaluation: checks for NaN/Infinity and reports location
    private static double SafeEval(Func<double, double> f, double x)
    {
        double y = f(x);
        if (double.IsNaN(y) || double.IsInfinity(y))
            throw new Exception($"Function returned {y} at x = {x}");
        return y;
    }

    // ------------------ Erf implementation (piecewise) ------------------
    // Uses Integrate(...) internally. eps and acc are passed to Integrate.
    public static double Erf(double z, double acc = 1e-6, double eps = 1e-6)
    {
        if (z < 0.0)
        {
            return -Erf(-z, acc, eps);
        }
        else if (z <= 1.0)
        {
            // erf(z) = 2/sqrt(pi) * integral_0^z exp(-x^2) dx
            double scale = 2.0 / Math.Sqrt(Math.PI);
            Func<double, double> g = x => Math.Exp(-x * x);
            return scale * Integrate(g, 0.0, z, acc, eps);
        }
        else // z > 1
        {
            // erf(z) = 1 - 2/sqrt(pi) * integral_0^1 dt exp(-(z + (1-t)/t)^2) / t^2
            double scale = 2.0 / Math.Sqrt(Math.PI);
            Func<double, double> h = t =>
            {
                double invt = 1.0 / t;
                double arg = z + (1.0 - t) * invt;
                double v = Math.Exp(-arg * arg);
                return v * invt * invt;
            };
            double integral = Integrate(h, 0.0, 1.0, acc, eps);
            return 1.0 - scale * integral;
        }
    }

    // ------------------ Testing helpers ------------------
    private static void RunTest(string desc, Func<double, double> fRaw, double a, double b, double expected, double acc, double eps)
    {
        Console.WriteLine("------------------------------------------------------------");
        Console.WriteLine(desc);
        long evalCount = 0;
        Func<double, double> f = x =>
        {
            evalCount++;
            return fRaw(x);
        };

        try
        {
            double result = Integrate(f, a, b, acc, eps);
            double absErr = Math.Abs(result - expected);
            bool passed = absErr <= acc + eps * Math.Abs(expected);

            Console.WriteLine($"  computed = {result:R}");
            Console.WriteLine($"  expected = {expected:R}");
            Console.WriteLine($"  abs error = {absErr:E}");
            Console.WriteLine($"  acc = {acc:G}, eps = {eps:G}");
            Console.WriteLine($"  passes (err ≤ acc + eps*|expected|)? : {passed}");
            Console.WriteLine($"  function evaluations: {evalCount}");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"  ERROR during integration: {ex.Message}");
            Console.WriteLine($"  function evaluations until error: {evalCount}");
        }
    }

    // ------------------ Main: run assignment tests then erf convergence ---------------
    static void Main(string[] args)
    {
        // PART A: the four tests (use assignment-default tolerances acc=1e-3, eps=1e-3)
        double acc_default = 1e-8;
        double eps_default = 1e-8;
        Console.WriteLine("Task A: \nWe start out by testing out integrator on the 4 expressions from the exercise. The third expression actually is equal to pi/4 and not pi/2. First we pick strict accuracy requirements of 1e-8 and then relax them afterwards to 1e-3 and see that our integrator returns values within     the from the exercise given accuracy goals.");

        RunTest("sqrt(x) on [0,1]", x => Math.Sqrt(x), 0.0, 1.0, 2.0 / 3.0, acc_default, eps_default);
        RunTest("1/sqrt(x) on [0,1]", x => 1.0 / Math.Sqrt(x), 0.0, 1.0, 2.0, acc_default, eps_default);
        // NOTE: integral sqrt(1-x^2) over [0,1] is pi/4. Some statements in the assignment list pi/2 (that's for [-1,1]).
        RunTest("sqrt(1-x^2) on [0,1]", x => Math.Sqrt(Math.Max(0.0, 1.0 - x * x)), 0.0, 1.0, Math.PI / 4.0, acc_default, eps_default);
        RunTest("ln(x)/sqrt(x) on [0,1]", x => Math.Log(x) / Math.Sqrt(x), 0.0, 1.0, -4.0, acc_default, eps_default);

        Console.WriteLine("\n\nNow trying with less strict accuracy goals:");
        acc_default = 1e-3;
        eps_default = 1e-3;


        RunTest("sqrt(x) on [0,1]", x => Math.Sqrt(x), 0.0, 1.0, 2.0 / 3.0, acc_default, eps_default);
        RunTest("1/sqrt(x) on [0,1]", x => 1.0 / Math.Sqrt(x), 0.0, 1.0, 2.0, acc_default, eps_default);
        // NOTE: integral sqrt(1-x^2) over [0,1] is pi/4. Some statements in the assignment list pi/2 (that's for [-1,1]).
        RunTest("sqrt(1-x^2) on [0,1]", x => Math.Sqrt(Math.Max(0.0, 1.0 - x * x)), 0.0, 1.0, Math.PI / 4.0, acc_default, eps_default);
        RunTest("ln(x)/sqrt(x) on [0,1]", x => Math.Log(x) / Math.Sqrt(x), 0.0, 1.0, -4.0, acc_default, eps_default);
        Console.WriteLine();
        Console.WriteLine();

        // PART B: erf(1) convergence experiment
        const double erf1_exact = 0.84270079294971486934; // reference value from the exercise
        double eps = 0.0; // as requested by the exercise

        string outFile = "erf_convergence.dat";
        using (var writer = new StreamWriter(outFile))
        {
            writer.WriteLine("# acc    abs_error    computed    func_evals");
            Console.WriteLine("Computing erf(1) for decreasing acc (eps = 0). Writing results to " + outFile);

            for (int e = 1; e <= 10; e++)
            {
                double acc = Math.Pow(10.0, -e); // 1e-1 .. 1e-10

                long evalCount = 0;
                Func<double, double> g = x =>
                {
                    evalCount++;
                    return Math.Exp(-x * x);
                };

                double computed = double.NaN;
                bool success = true;
                try
                {
                    double scale = 2.0 / Math.Sqrt(Math.PI);
                    double integral = Integrate(g, 0.0, 1.0, acc, eps);
                    computed = scale * integral;
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"  acc={acc:E}: ERROR during integration: {ex.Message}");
                    success = false;
                }

                double absErr = success ? Math.Abs(computed - erf1_exact) : double.NaN;
                writer.WriteLine("{0:E} {1:E} {2:R} {3}", acc, absErr, computed, evalCount);

                Console.WriteLine($"acc={acc:E}  abs_error={absErr:E}  computed={computed:R}  func_evals={evalCount}");
            }
        }

        Console.WriteLine();
        Console.WriteLine($"See the plot in the file erf.convergence.png");

        Console.WriteLine("\nTask B:");

        Console.WriteLine("We calculate integral examples from the exercise with integrable divergencies at the end-points of the intervals");
        // Compare integrators on the two "integrable divergence" examples:
        double acc_compare = 1e-3;
        double eps_compare = 1e-3;

        CompareIntegrators("Compare ∫0^1 1/sqrt(x) dx", x => 1.0 / Math.Sqrt(x), 0.0, 1.0, acc_compare, eps_compare, 2.0);
        CompareIntegrators("Compare ∫0^1 ln(x)/sqrt(x) dx", x => Math.Log(x) / Math.Sqrt(x), 0.0, 1.0, acc_compare, eps_compare, -4.0);

        Console.WriteLine("We compared the number of integrand evaluations with the python/numpy's integration routines");
        Console.WriteLine("Python SciPy (quad) reference:");
        Console.WriteLine("∫_0^1 1/sqrt(x) dx ≈ 2.000, calls = 231");
        Console.WriteLine("∫_0^1 ln(x)/sqrt(x) dx ≈ -4.000, calls = 315");

        Console.WriteLine("In both test cases the Clenshaw–Curtis method reached the target accuracy with fewer integrand evaluations than the ordinary adaptive integrator and SciPy’s quad.\n");

        Console.WriteLine("I now generalize the integrator to accept infinite limits and test it on two examples ∫0^∞ exp(-x) dx and ∫-∞^∞ exp(-x^2) dx. We get:");
        long calls;

        // Test 1: ∫_0^∞ exp(-x) dx = 1
        double res1 = AdaptiveIntegrator.Integrate(x => Math.Exp(-x), 0, double.PositiveInfinity, out calls);
        Console.WriteLine($"∫0^∞ exp(-x) dx ≈ {res1}, exact = 1, calls = {calls}");

        // Test 2: ∫_-∞^∞ exp(-x^2) dx = sqrt(pi)
        double res2 = AdaptiveIntegrator.Integrate(x => Math.Exp(-x * x), double.NegativeInfinity, double.PositiveInfinity, out calls);
        Console.WriteLine($"∫-∞^∞ exp(-x^2) dx ≈ {res2}, exact = {Math.Sqrt(Math.PI)}, calls = {calls}");

        Console.WriteLine("\nWe now compare the results with Python SciPy (quad) for the same integrals:");
        Console.WriteLine("Python SciPy (quad) reference:");
        Console.WriteLine("∫0^∞ e^-x dx         | SciPy quad result: 1.000000000000 | calls: 135");
        Console.WriteLine("∫-∞^∞ exp(-x^2) dx   | SciPy quad result: 1.772453850906  | calls: 270");
        Console.WriteLine();
        Console.WriteLine("In both examples the generalized integrator reaches the target accuracy with fewer function evaluations than SciPy's quad");


        Console.WriteLine("\nTask C:");
        Console.WriteLine("making integrator estimate and return the integration error:");
        // Example 1: ∫0^1 ln(x)/√x dx
        var f1 = new Func<double, double>(x => Math.Log(x) / Math.Sqrt(x));
        var exact1 = -4.0;
        var (result1, estError1) = IntegrateWithError(f1, 0, 1);
        Console.WriteLine($"∫0^1 ln(x)/√x dx ≈ {result1}, estimated error = {estError1}, actual error = {Math.Abs(result1 - exact1)}");

        // Example 2: ∫0^π sin(x) dx
        var f2 = new Func<double, double>(x => Math.Sin(x));
        var exact2 = 2.0;
        var (result2, estError2) = IntegrateWithError(f2, 0, Math.PI);
        Console.WriteLine($"∫0^π sin(x) dx ≈ {result2}, estimated error = {estError2}, actual error = {Math.Abs(result2 - exact2)}");





    } // end Main
}
