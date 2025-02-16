using static System.Console;
using static System.Math;
using sfuns;
class math{

static int Main(){
    Write("\nTask 1:\n");
    double sqrt2 = Sqrt(2.0);
    WriteLine($"Sqrt(2) = {sqrt2}");
    //double sqrt2=System.Math.Sqrt(2.0);
    //WriteLine($"Sqrt(2) = {sqrt2}");
    double sqrt22 = Pow(2,1.0/5);
    WriteLine($"5th root of 2 = {sqrt22}");
    double jj = Pow(E,PI);
    WriteLine($"e^pi = {jj}");
    jj = Pow(PI,E);
    WriteLine($"Pi^e = {jj}\n");
    WriteLine("Task 2:");
    for (int i = 1; i <= 10; i++)
        {
            double gamma = Sfuns.fgamma(i);
            WriteLine($"*Γ({i}) = {gamma}");
        }
    WriteLine("the actual factorial results:");
    for (int i = 1; i < 11; i++)
    {
        
        double Gamme_n = 1;
        for (int j = 1; j < i; j++)
        {
            Gamme_n *= j;
        }
        WriteLine($"Γ({i}) = {Gamme_n}");

    }
    WriteLine("They are close enough\n");
    
    WriteLine("Task 3:");
    WriteLine("The result of lngamma for n = 1,2,...,10");
    for (int i = 1; i <= 10; i++)
        {
            double gamma = Sfuns.lngamma(i);
            WriteLine($"*Γ({i}) = {gamma}");
        }
	return 0;
    
	}
}