using static System.Math;

namespace sfuns
{
    public static class Sfuns
    {
        public static double fgamma(double x)
        {
            if (x < 0) return PI / Sin(PI * x) / fgamma(1 - x); // Euler's reflection formula
            if (x < 9) return fgamma(x + 1) / x; // Recurrence relation
            double lnfgamma = x * Log(x + 1 / (12 * x - 1 / x / 10)) - x + Log(2 * PI / x) / 2;
            return Exp(lnfgamma);
        }
        public static double lngamma(double x)
        {
            if (x < 0) return double.NaN; // Euler's reflection formula
            if (x < 9) return lngamma(x + 1) - Log(x); // Recurrence relation
            return x * Log(x + 1 / (12 * x - 1 / x / 10)) - x + Log(2 * PI / x) / 2;
        }
    }
}
