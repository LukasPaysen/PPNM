using static System.Console;

class epsilon{

static int Main(){
    Write("\nTask 1:\n");
    int i=1; 
    while(i+1>i)
    {
        i++;
    } 
    Write("my max int = {0}\n",i);
    Write("int.MaxValue = {0}\n",int.MaxValue);
    int j = 1;
    while(j-1<j)
    {
        j--;
    }
    Write("my min int = {0}\n",j);
    Write("int.MinValue = {0}\n\n", int.MinValue);
    Write("Task 2: \n");

    double x=1; 
    while(1+x!=1)
    {x/=2;} 
    
    x *= 2;
    WriteLine("machine Epsilon double = {0}", x); 
    float y=1F; 
    while((float)(1F+y) != 1F)
    {
        y/=2F;
    } 
    y*=2F;
    WriteLine("machine epsilon float = {0}\n", y);
    WriteLine("doble epsilon should be about {0}", System.Math.Pow(2,-52));
    WriteLine("float epsilon should be about {0}", System.Math.Pow(2,-23));
    WriteLine("Checks out\n");
    WriteLine("Task 3:");
    double epsilon=System.Math.Pow(2,-52);
    double tiny=epsilon/2;
    double a=1+tiny+tiny;
    double b=tiny+tiny+1;
    WriteLine("{0}, {1}", a, b);
    Write($"a==b ? {a==b}\n");
    Write($"a>1  ? {a>1}\n");
    Write($"b>1  ? {b>1}\n");
    Write("This happens because c# first adds the two first elements together, which for b is tiny + tiny = epsilon and then adds 1 => b=1+epsilon, whereas a+tiny =1 and then 1+ tiny=a=1\n \n");

    WriteLine("Task 4:");
    double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
    double d2 = 8*0.1;
    WriteLine($"d1={d1:e15}");
    WriteLine($"d2={d2:e15}");
    WriteLine($"d1==d2 ? => {d1==d2}"); 
    WriteLine("are d1 and d2 equal? We check using approx, approx(a,b) = {0}",approx(d1,d2));

    
	return 0;
    
	}public static bool approx(double a, double b, double acc=1e-9, double eps=1e-9)
    {
        if(System.Math.Abs(b-a) <= acc) return true;
        if(System.Math.Abs(b-a) <= System.Math.Max(System.Math.Abs(a),System.Math.Abs(b))*eps) return true;
        return false;
    }
}