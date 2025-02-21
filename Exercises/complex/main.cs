using System;
using System.Numerics;

class Program
{
    // Returns true if a and b are approximately equal.
    static bool Approx(double a, double b, double acc = 1e-9, double eps = 1e-9)
    {
        if (Math.Abs(a - b) < acc) return true;
        if (Math.Abs(a - b) < (Math.Abs(a) + Math.Abs(b)) * eps) return true;
        return false;
    }
    
    // Compares two (real, imaginary) pairs representing complex numbers.
    static bool AreTuplesEqual(Tuple<double, double> t1, Tuple<double, double> t2, double acc = 1e-9, double eps = 1e-9)
    {
        return Approx(t1.Item1, t2.Item1, acc, eps) && Approx(t1.Item2, t2.Item2, acc, eps);
    }
    
    // Helper to compare two Complex numbers.
    static bool AreComplexEqual(Complex c1, Complex c2, double acc = 1e-9, double eps = 1e-9)
    {
        return AreTuplesEqual(new Tuple<double, double>(c1.Real, c1.Imaginary), 
                              new Tuple<double, double>(c2.Real, c2.Imaginary), 
                              acc, eps);
    }
    
    // Compares a group of computed values and prints them.
    // It prints each label and value, then checks that all values are pairwise equal.
    static void CompareGroup(string description, Tuple<string, Complex>[] items)
    {
        bool allEqual = true;
        for (int i = 1; i < items.Length; i++)
        {
            if (!AreComplexEqual(items[0].Item2, items[i].Item2))
            {
                allEqual = false;
                break;
            }
        }
        
        Console.WriteLine($"{description}:");
        for (int i = 0; i < items.Length; i++)
        {
            Console.Write($"{items[i].Item1} = {items[i].Item2}");
            if (i < items.Length - 1)
                Console.Write(", ");
        }
        Console.WriteLine($". Equal? {allEqual}\n");
    }
    
    static void Main()
    {
        Complex i = Complex.ImaginaryOne;
        
        // Group 1: √-1 comparisons.
        Complex sqrtMinusOne = Complex.Sqrt(-1);
        Console.WriteLine("First we calculate:");
        Console.WriteLine("√-1 = " + Complex.Sqrt(-1));
        Console.WriteLine("√i = " + Complex.Sqrt(i));
        Console.WriteLine("e^i = " + Complex.Exp(i));
        Console.WriteLine("e^(iπ) = " + Complex.Exp(i * Math.PI));
        Console.WriteLine("i^i = " + Complex.Pow(i, i));
        Console.WriteLine("ln(i) = " + Complex.Log(i));
        Console.WriteLine("sin(iπ) = " + Complex.Sin(i * Math.PI) + "\n");
        Console.WriteLine("Now we compare if they are equal using our approx method. \n A method AreTuplesEqual is used to take two tuples and check if they are equal by checking each element through our Approx() method");
        // Define i.
        
        // Compare computed √-1 with i.
        CompareGroup("√-1 = i?", new Tuple<string, Complex>[] {
            new Tuple<string, Complex>("√-1", sqrtMinusOne),
            new Tuple<string, Complex>("i", i)
        });
        // Also compare computed √-1 with -i.
        CompareGroup("√-1 = -i?", new Tuple<string, Complex>[] {
            new Tuple<string, Complex>("√-1", sqrtMinusOne),
            new Tuple<string, Complex>("-i", -i)
        });
        
        // Group 2: ln(i) comparisons.
        // ln(i) should be iπ/2 and ln(e^(iπ/2)) should also be iπ/2.
        Complex ln_i = Complex.Log(i);
        Complex ln_exp = Complex.Log(Complex.Exp(new Complex(0, Math.PI / 2)));
        Complex expectedLn = new Complex(0, Math.PI / 2);
        CompareGroup("ln(i) = ln(e^(iπ/2)) = iπ/2?", new Tuple<string, Complex>[] {
            new Tuple<string, Complex>("ln(i)", ln_i),
            new Tuple<string, Complex>("ln(e^(iπ/2))", ln_exp),
            new Tuple<string, Complex>("iπ/2", expectedLn)
        });
        
        // Group 3: √i comparisons.
        // √i should be e^(iπ/4) = 1/√2 + i/√2.
        Complex sqrt_i = Complex.Sqrt(i);
        Complex expectedSqrtI = new Complex(Math.Cos(Math.PI / 4), Math.Sin(Math.PI / 4));
        CompareGroup("√i = e^(iπ/4) = 1/√2 + i/√2", new Tuple<string, Complex>[] {
            new Tuple<string, Complex>("√i", sqrt_i),
            new Tuple<string, Complex>("e^(iπ/4)", expectedSqrtI)
        });
        
        // Group 4: i^i comparisons.
        // i^i should equal e^(i ln(i)) = e^(-π/2) (a real number).
        Complex ii = Complex.Pow(i, i);
        Complex expectedII = new Complex(Math.Exp(-Math.PI / 2), 0);
        CompareGroup("i^i = e^(i ln(i)) = e^(-π/2)?", new Tuple<string, Complex>[] {
            new Tuple<string, Complex>("i^i", ii),
            new Tuple<string, Complex>("e^(i ln(i))", expectedII)
        });
    }
}
