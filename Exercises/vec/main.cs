using static System.Console;
using static System.Math;

public class vec{
	public double x,y,z; /* the three components of a vector */
	public vec(){ x=y=z=0; } // default constructor
	public vec(double x,double y,double z){ // parametrized constructor
		this.x=x; this.y=y; this.z=z;
		}
	/* ... more stuff goes here ... */
    public void print(string s=""){Write(s);WriteLine($"({x}, {y}, {z})");}
    public static vec operator*(vec v, double c){return new vec(c*v.x,c*v.y,c*v.z);}
    public static vec operator*(double c, vec v){return v*c;}
    public static vec operator/(vec u, double c){return new vec(u.x/c,u.y/c,u.z/c);}
    public static vec operator+(vec u, vec v){return new vec(u.x+v.x,u.y+v.y,u.z+v.z);}
    public static vec operator-(vec u){return new vec(-u.x,-u.y,-u.z);}
    public static vec operator-(vec u, vec v){return new vec(u.x-v.x,u.y-v.y,u.z-v.z);}

    public double dot(vec other) /* to be called as u.dot(v) */
	{return this.x*other.x+this.y*other.y+this.z*other.z;}
    public static double dot(vec v,vec w) /* to be called as vec.dot(u,v) */
	{return v.x*w.x+v.y*w.y+v.z*w.z;}

    static bool approx(double a,double b,double acc=1e-9,double eps=1e-9){
	if(Abs(a-b)<acc)return true;
	if(Abs(a-b)<(Abs(a)+Abs(b))*eps)return true;
	return false;
	}

    public override string ToString(){ return $"{x} {y} {z}"; }
    public bool approx(vec other){
	if(!approx(this.x,other.x))return false;
	if(!approx(this.y,other.y))return false;
	if(!approx(this.z,other.z))return false;
	return true;
	}
    public static bool approx(vec u, vec v) => u.approx(v);

    static int Main()
    {
        WriteLine("Task 1: \n");
        vec u = new vec(3,3,3);
        vec v = new vec(2,2,2);
        u.print("vector u = ");
        v.print("vector v = ");

        WriteLine("\nTask 2:");
        vec x = u/3;
        x.print("testing /: dividing u by 3: u/3 = ");
        x = u +v;
        x.print("testing +: u + v = ");
        x = u-v;
        x.print("testing -: u - v = ");

        WriteLine("\nTask 3:");

        
        double k = u.dot(v);
        WriteLine("u.dot(v) = {0}", k);
        k = vec.dot(u,v);
        WriteLine("vec.dot(u,v) = {0}", k);
        
        WriteLine($"Are u and v equal according to our approx method? u.approx(v) = {u.approx(v)}");
        u = new vec(2,2,2);
        WriteLine($"What if u = (2,2,2)? u.approx(v) = {u.approx(v)}");
        
        WriteLine("\nTask 4: ");
        WriteLine("Testing our overwritten ToString() method: u.ToString() = {0}", u.ToString());

        
        return 0;  

    }
}



