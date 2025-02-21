using System;
using static System.Console;
using static System.Math;

class Program {
    static void Main() {
        char[] split_delimiters = { ' ', '\t', '\n' };
        var split_options = StringSplitOptions.RemoveEmptyEntries;

        for (string line = ReadLine(); line != null; line = ReadLine()) {
            var numbers = line.Split(split_delimiters, split_options);
            foreach (var number in numbers) {
                double x = double.Parse(number);
                Error.WriteLine($"{x} {Sin(x)} {Cos(x)}");
            }
        }
    }
}
