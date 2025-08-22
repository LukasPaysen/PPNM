Practical Programning and Numerical Methods: Exam 2025

Least-Squares Signal Declipping

1.
My implementation is in main.cs (single file). It contains small linear-algebra helpers (vector, matrix), a QR solver (classical Gram–Schmidt, tall matrices), and a Declipping class with:

 - BuildThirdDerivativeD(N): builds the finite-difference matrix D for the 3rd derivative.
Interior rows use the 5-point stencil −½,1,0,−1,½; the first/last rows use forward/backward −1,3,−3,1 style stencils.

 - Declip(vector y, double y_min, double y_max, bool enforceConsistency=true): core routine.

The method follows the note: write the reconstructed signal as

x = ỹ + M z


where ỹ is y with clipped entries set to 0, M inserts the unknowns z back into those positions.
We minimize ‖D(ỹ + Mz)‖² which gives the least-squares system

A z = −D ỹ,   with  A = D M .


I form A by taking the columns of D at the clipped indices and solve with my own QR factorization.

 How to use
Call:
var xDeclipped = Declipping.Declip(y, y_min, y_max);
y is the clipped vector; the function returns the reconstructed vector x.

2. (Sine test)
Generated a sine with 5 cycles on 0,1 and clipped at ±0.8 (N=1000).

Data written to declipping.dat (columns: i, t, x_true, y_clipped, x_declipped, clipped_flag).

Plot saved as declipping.png.
Result: the declipped curve is basically on top of the true curve; the printed summary shows max |error| < 2e-4.

3. (Harder signal)
A “made-up” signal: chirp + harmonic + slow baseline wobble + a short burst, then clipped at ±0.7 (N=2000).

Data in declipping_complex.dat; plot declipping_complex.png.
Outcome: still good overall—largest mismatch is across a long fully-clipped block (~0.06–0.09 s) where the method smoothly bridges and overshoots a bit.

I tested both tasks from the exam (simple synthetic sine and a more complicated example) and included the output files and plots.
