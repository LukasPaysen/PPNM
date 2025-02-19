using System;
using System.IO;
using static System.Console;
using static System.Math;

class Program {
    static int Main(string[] args) {
        string infile = null, outfile = null;

        // Parse command-line arguments
        foreach (string arg in args) {
            string[] words = arg.Split(':');
            if (words[0] == "-input") infile = words[1];
            if (words[0] == "-output") outfile = words[1];
        }

        // Check if both input and output files are specified
        if (infile == null || outfile == null) {
            Error.WriteLine("wrong filename argument");
            return 1; // Return an error code
        }

        try {
            // Open input and output streams
            StreamReader instream = new StreamReader(infile);
            StreamWriter outstream = new StreamWriter(outfile, append: false);

            // Read numbers from input file and write their sine and cosine
            string line;
            while ((line = instream.ReadLine()) != null) {
                double x = double.Parse(line);  // Parse the number from the line
                outstream.WriteLine($"{x} {Sin(x)} {Cos(x)}");  // Write number, sine, and cosine to output file
            }

            instream.Close();
            outstream.Close();
        }
        catch (Exception ex) {
            Error.WriteLine($"Error: {ex.Message}");
            return 1;
        }

        return 0; // Success
    }
}
