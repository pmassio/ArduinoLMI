# ArduinoLMI
Linear Matrix Inequality (LMI) solver for Arduino. Allows the resolution of a set of simple standard LMI problems such as robust pole placement, robust stabilisation, gain-scheduled control.
## References
If you use this code for research purposes, please cite the article  [*Hiba Houmsi, Paolo Massioni, Federico Bribiesca Argomedo, Romain Delpoux. "Embedded controller optimization for efficient electric motor drive," 2023 IEEE Vehicle Power and Propulsion (VPPC 2023), Oct 2023, Milan, Italy*](https://ec-lyon.hal.science/hal-04176290/). The article contains a description of the embedded implementation of the interior-point method solver used in this code.

For a general introduction on Linaear Matrix inequalities (LMI) problems, refer to 
[*Stephen Boyd, Laurent El Ghaoui, E. Feron, 
 and V. Balakrishnan, "Linear Matrix Inequalities in System and Control Theory",
 Volume 15 of Studies in Applied Mathematics, Society for Industrial and Applied 
 Mathematics (SIAM), 1994*](https://web.stanford.edu/~boyd/lmibook/lmibook.pdf),
which also provides a description of the interior point algorithm.

For the regional pole placement method, refer to 
 [*M. Chilali, P. Gahinet and P. Apkarian, "Robust pole 
 placement in LMI regions," IEEE Transactions on Automatic Control, vol. 44, no. 12,
 pp. 2257-2270, Dec. 1999*](https://doi.org/10.1109/9.811208).

## Requirements
Requires [ArduinoEigen](https://github.com/hideakitai/ArduinoEigen), so it will work only on architectures supporting c++ standard libraries (STM and SAM). It might work on other architectures if you manage to install such libraries somehow, with a bit more effort. 

## Usage
The library provides an embedded implementation of controller synthesis problems requiring the solution of a Linear Matrix Inequality (LMI) feasibility problem. In general this can be made on a laptop computer, and the resulting gain can be just copied on the embedded controller; this library is aimed at applications where the controller needs to be recomputed on board, such as adaptive control, autonomous systems, etc. Please take into account that the solver will take a significant time to deliver a solution; I have tested it on an MKR1010 and on a GIGA R1; for the first one expect a few seconds or few tens of seconds for systems up to order four, whereas for the GIGA R1 it runs a lot faster. 

Do not use this code for high-order system.

The functions return $0$ matrices if the problem is unfeasible. There are no checks on the input matrices size, make sure the $A$ matrices are square, and that the $B$ matrices have as many rows as the $A$ matrices.

Please see the example.ino for an example of application code.

## User functions

### Nominal decay rate stabilisation
```cpp
Eigen::MatrixXf DecayRate(Eigen::MatrixXf& A, Eigen::MatrixXf& B, float amin)  
```
Returns a MatriXf object containing a matrix $K$ of gains, such that the matrix 
$A+BK$
is Hurwitz and has all the eigenvalues $\lambda$ fulfilling the following constraint:
<br> $\Re(\lambda) \leqslant -\texttt{amin}$. 
<br>
This is equivalent to imposing a minimum (positive) decay rate  $\texttt{amin}$.


### Nominal regional pole placement
```cpp
Eigen::MatrixXf RegionalPolePlacement(Eigen::MatrixXf& A, Eigen::MatrixXf& B, float amax, float amin, float beta)  
```
Returns a MatriXf object containing a matrix $K$ of gains, such that the matrix 
$A+BK$
is Hurwitz and has all the eigenvalues $\lambda$ fulfilling the following constraints:
<br> $-\texttt{amax} \leqslant \Re(\lambda) \leqslant -\texttt{amin}$, 
<br> $|\Im(\lambda)|\leqslant \texttt{beta} |\Re(\lambda)|$.
<br>
This is equivalent to imposing a minimum (positive) decay rate  $\texttt{amin}$, a maximum decay rate  $\texttt{amax}$, and a minimum damping rate.

### Robust decay rate stabilisation (with 1 parameter)
```cpp
Eigen::MatrixXf DeacayRate1Param(Eigen::MatrixXf& A0, Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amin,  float p1max, float p1min);    
```
Returns a MatriXf object containing a matrix $K$ of gains, such that the matrix 
$A_0+p_1 A_1+ BK$
is Hurwitz and has all the eigenvalues $\lambda$ fulfilling the following constraint:
<br> $\Re(\lambda) \leqslant -\texttt{amin}$, 
<br>
for all values of the parameter $p_1$ in the interval $[\texttt{p1max}, \texttt{p1min}]$.
This is equivalent to imposing a minimum (positive) decay rate  $\texttt{amin}$. This property is kept for any arbitrarily fast change of the value of the parameter (within the interval).

### Robust regional pole placement (with 1 parameter)
```cpp
Eigen::MatrixXf RegionalPolePlacementRobust1Param(Eigen::MatrixXf& A0, Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amax, float amin, float beta, float p1max, float p1min);    
```
Returns a MatriXf object containing a matrix $K$ of gains, such that the matrix 
$A_0+p_1 A_1+ BK$
is Hurwitz and has all the eigenvalues $\lambda$ fulfilling the following constraints:
<br> $-\texttt{amax} \leqslant \Re(\lambda) \leqslant -\texttt{amin}$, 
<br> $|\Im(\lambda)|\leqslant \texttt{beta} |\Re(\lambda)|$,
<br>
for all values of the parameter $p_1$ in the interval $[\texttt{p1max}, \texttt{p1min}]$.
This is equivalent to imposing a minimum (positive) decay rate  $\texttt{amin}$, a maximum decay rate  $\texttt{amax}$, and a minimum damping rate. These properties are kept for any arbitrarily fast change of the value of the parameter (within the interval).

### Robust decay rate stabilisation under norm-constrained bound
```cpp
Eigen::MatrixXf DecayRateRobustNormBound(Eigen::MatrixXf& A0, Eigen::MatrixXf& B, Eigen::MatrixXf& E, Eigen::MatrixXf& F, float amin);  
```
Returns a MatriXf object containing a matrix $K$ of gains, such that the matrix 
$A_0+E \Delta F+ BK$
is Hurwitz and has all the eigenvalues $\lambda$ fulfilling the following constraint:
<br> $\Re(\lambda) \leqslant -\texttt{amin}$, 
<br>
for all values of the $d \times d$ matrix parameter $\Delta$ with maximum singular value smaller or equal to 1.
This is equivalent to imposing a minimum (positive) decay rate  $\texttt{amin}$. This property is kept for any arbitrarily fast change of the value of the parameter (within the norm constraint). Notice that $E \in \matbb{R}^{n \times dn}$ and $F \in \matbb{R}^{d \times n}$ (with $A \in \matbb{R}^{n \times n}$); make sure that these parameters are entered correctly because there is no internal check on their dimensions or consistency.
  

### Gain-scheduled decay rate stabilisation (with 1 parameter)
```cpp
Eigen::MatrixXf** DecayRateGainScheduling1Param(Eigen::MatrixXf& A0,  Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amin, float pmax, float pmin);
```
Returns a pointer to a two-dimensional array of MatriXf pointers object containing a matrix $K$ of gains, such that the matrix 
$A_0+p_1 A_1+ B(\texttt{*K[0]}+p_1\texttt{*K[1]})$
is Hurwitz and has all the eigenvalues $\lambda$ fulfilling the following constraint:
<br> $\Re(\lambda) \leqslant -\texttt{amin}$, <br>
for all values of the parameter $p_1$ in the interval $[\texttt{p1max}, \texttt{p1min}]$.
This is equivalent to imposing a minimum (positive) decay rate  $\texttt{amin}$. This property is kept for any arbitrarily fast change of the value of the parameter (within the interval). Remember to delete  \texttt{Kp[0]}
  and \texttt{Kp[1]} to free the memory before launching this function for a second time.


### Gain-scheduled regional pole placement (with 1 parameter)
```cpp
Eigen::MatrixXf** RegionalPolePlacementGainScheduling1Param(Eigen::MatrixXf& A0,  Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amax, float amin, float beta, float pmax, float pmin);
```
Returns a pointer to a two-dimensional array of MatriXf pointers object containing a matrix $K$ of gains, such that the matrix 
$A_0+p_1 A_1+ B(\texttt{*K[0]}+p_1\texttt{*K[1]})$
is Hurwitz and has all the eigenvalues $\lambda$ fulfilling the following constraints:
<br> $-\texttt{amax} \leqslant \Re(\lambda) \leqslant -\texttt{amin}$, 
<br> $|\Im(\lambda)|\leqslant \texttt{beta} |\Re(\lambda)|$,
<br>
for all values of the parameter $p_1$ in the interval $[\texttt{p1max}, \texttt{p1min}]$.
This is equivalent to imposing a minimum (positive) decay rate  $\texttt{amin}$, a maximum decay rate  $\texttt{amax}$, and a minimum damping rate. These properties are kept for any arbitrarily fast change of the value of the parameter (within the interval). Remember to delete  \texttt{Kp[0]}
  and \texttt{Kp[1]} to free the memory before launching this function for a second time.




# Licence 
GNU General Public License 3.0


# Autor
[Paolo Massioni](https://sites.google.com/site/pmassio/home)

