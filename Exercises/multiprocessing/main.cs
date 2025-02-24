using System;
using System.Threading;
using System.Threading.Tasks;
using System.Linq;
public class Data {
    public int a, b;
    public double sum;
}

class Program {
    public static void Harm(object obj) {
        var arg = (Data)obj;
        arg.sum = 0;
        for (int i = arg.a; i < arg.b; i++) {
            arg.sum += 1.0 / i;
        }
    }

    static void Main(string[] args) {
        // Default values
        int nthreads = 1, nterms = (int)1e8;

        // Parse command-line arguments
        foreach (var arg in args) {
            var words = arg.Split(':');
            if (words[0] == "-threads") nthreads = int.Parse(words[1]);
            if (words[0] == "-terms") nterms = (int)float.Parse(words[1]);
        }
        if (nthreads == 1) {System.Console.WriteLine("\nTask 1:");}
        if (nthreads < 5) 
        {
        // Prepare data objects
        Data[] parameters = new Data[nthreads];
        for (int i = 0; i < nthreads; i++) {
            parameters[i] = new Data();
            parameters[i].a = 1 + nterms / nthreads * i;
            parameters[i].b = 1 + nterms / nthreads * (i + 1);
        }
        parameters[parameters.Length - 1].b = nterms + 1; // Adjust endpoint

        // Create and start threads
        Thread[] threads = new Thread[nthreads];
        for (int i = 0; i < nthreads; i++) {
            threads[i] = new Thread(Harm);
            threads[i].Start(parameters[i]);
        }

        // Join threads (wait for them to finish)
        foreach (var thread in threads) thread.Join();

        // Calculate total sum
        double total = 0;
        foreach (var param in parameters) total += param.sum;


        // Print result
        Console.WriteLine($"Harmonic sum ({nterms} terms) with {nthreads} threads: {total}");
        }













        if(nthreads == 5) 
        {
        System.Console.WriteLine("\nTask 2:");

        double serialSum = 0;
        for (int i = 1; i <= nterms; i++)
        {
            serialSum += 1.0 / i;
        }
        Console.WriteLine($"Serial sum ({nterms} terms): {serialSum}");
        }
        if(nthreads == 6) {
        double parallelSum = 0;

        // Parallel.For loop without synchronization (race condition)
        Parallel.For(1, nterms + 1, (i) =>
        {
            parallelSum += 1.0 / i; // This causes a race condition
        });

        var parallelEndTime = DateTime.Now;
        Console.WriteLine($"Parallel sum ({nterms}): {parallelSum}");
        
        
        }
        else if(nthreads == 7)
        {
            Console.WriteLine("When you use too many threads for a small task, the time taken to manage these threads outweighs the parallelization benefit.");
        Console.WriteLine("Incorrect result is due to race conditions.");
        var sum = new System.Threading.ThreadLocal<double>(() => 0, trackAllValues: true);

        // Run the parallel for-loop
        Parallel.For(1, nterms + 1, (i) =>
        {
            sum.Value += 1.0 / i;
        });

        // Sum all thread-local values
        double totalSum = sum.Values.Sum();
        Console.WriteLine($"Total sum using ThreadLocal: {totalSum}");
        }
        else if (nthreads == 8)
        {
            Console.WriteLine("The use of ThreadLocal<T> gives the correct result but is slower due to the overhead of managing multiple threads, context switching, and synchronization");
        }
        }
        
}



