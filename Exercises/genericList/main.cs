
using System;
using static System.Console;
using static System.Math;


class Program {
    static void Main() {
        Write("Task 1: \n")
        var list = new genlist<double[]>(); 
        char[] delimiters = {' ','\t'};
        var options = StringSplitOptions.RemoveEmptyEntries;
        for(string line = ReadLine(); line!=null; line = ReadLine()){
            var words = line.Split(delimiters,options);
            int n = words.Length;
            var numbers = new double[n];
            for(int i=0;i<n;i++) numbers[i] = double.Parse(words[i]);
            list.add(numbers);
                }
        for(int i=0;i<list.size;i++){
            var numbers = list[i];
            foreach(var number in numbers)Write($"{number : 0.00e+00;-0.00e+00} ");
            WriteLine();
                }

    }
}
public class genlist<T>{
	public T[] data;
	public int size => data.Length; // property
	public T this[int i] => data[i]; // indexer
	public genlist(){ data = new T[0]; }
	public void add(T item){ /* add item to the list */
		T[] newdata = new T[size+1];
		System.Array.Copy(data,newdata,size);
		newdata[size]=item;
		data=newdata;
	}
}