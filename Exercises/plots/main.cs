using static System.Console;
using static System.Math;
using System;


using System.Collections.Generic;
class program{
    static void Main()
    {
        List<double> secondColumnValues = new List<double>
        {
            0.022564575,
            0.045111106,
            0.067621594,
            0.090078126,
            0.112462916,
            0.222702589,
            0.328626759,
            0.428392355,
            0.520499878,
            0.603856091,
            0.677801194,
            0.742100965,
            0.796908212,
            0.842700793,
            0.880205070,
            0.910313978,
            0.934007945,
            0.952285120,
            0.966105146,
            0.976348383,
            0.983790459,
            0.989090502,
            0.992790429,
            0.995322265
        };
        var firstColumn = new List<float> { 0.02f, 0.04f, 0.06f, 0.08f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2f };

        
        for (int i = 0; i < 100; i++)
        {
            float j = (float)i/50;
            float k = ((float)i +10)/10;
            if (i<24)
            {
            //float j = (float)program.secondColumnValues[i];
            float result = 1;
            float p = (float)((i+2.5)/2.5);
            for (float n = p; n > 1; n--)
            {
                result *= n;

            }
            float h = (float)Log(result);
            WriteLine(j.ToString() + " " + k.ToString() + " " + erf(j).ToString() + " " + sgamma(k+1).ToString() + " " + lngamma(k+1) + " " + p.ToString() + " " + result.ToString() + " " + h.ToString() + " " + i.ToString() + " " + secondColumnValues[i].ToString() + " " + firstColumn[i].ToString());
            }

        

            WriteLine(j.ToString() + " " + k.ToString() + " " + erf(j).ToString() + " " + sgamma(k+1).ToString() + " " + lngamma(k+1));
        }/*

        
        for (int i = 0; i < 24; i++)
        {
            //float j = (float)program.secondColumnValues[i];
            float result = 1;
            float j = (float)((i+2.5)/2.5);
            for (float n = j; n > 1; n--)
            {
                result *= n;

            }
            float h = (float)Log(result);
            WriteLine(j.ToString() + " " + result.ToString() + " " + h.ToString() + " " + i.ToString() + " " + secondColumnValues[i].ToString());

        }*/
    }

    

static double erf(double x){
/// single precision error function (Abramowitz and Stegun, from Wikipedia)
if(x<0) return -erf(-x);
double[] a={0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
double t=1/(1+0.3275911*x);
double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
return 1-sum*Exp(-x*x);
} 
public static double sgamma(double x){
if(x<0)return PI/Sin(PI*x)/sgamma(1-x);
if(x<9)return sgamma(x+1)/x;
double lnsgamma=Log(2*PI)/2+(x-0.5)*Log(x)-x
    +(1.0/12)/x-(1.0/360)/(x*x*x)+(1.0/1260)/(x*x*x*x*x);
return Exp(lnsgamma);
}
static double lngamma(double x){
if(x<=0) throw new ArgumentException("lngamma: x<=0");
if(x<9) return lngamma(x+1)-Log(x);
return x*Log(x+1/(12*x-1/x/10))-x+Log(2*PI/x)/2;
}



}