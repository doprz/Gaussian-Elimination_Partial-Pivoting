# M 340L (Matrices and Matrix Calculations) Coding Project
Instructor: Dr. Axel G. R. Turnquist

Note: All code used in the project must be attached to the project file.

## Computational Problem
Code up row reduction for an m x n matrix. Test with an 11 x 10 matrix with random (uniformly distributed between 0 and 1) entries. Print the result of row reduction to echelon form.

*Commentary: For this project you do not have to check for zeros. The example you will compute will have numbers chosen from the uniform distribution between 0 and 1, so with zero probability will you have a zero anywhere in the matrix. However, your code should be written to row reduce a m x n matrix. You need to determine how you choose the pivots. That’s up to you.*

## Theoretical Problem
Explicitly count the floating point operations FLOPs required in turning a m×n matrix into echelon form. *Hint:* How to count flops. Each multiplication, division, subtraction, and addition counts as one floating point. Usually this number in the end is reported using order notation. That is, if we have $\mathcal{O}(n)$ (read: “order n”) it means in a nutshell that the process scales no worse than some constant times n as the number n increases to infinity.

*Commentary: There are many resources for finding out how to count FLOPS online. Generally, Gaussian elimination (row reduction) for an n x n matrix is $\mathcal{O}(n^3)$ (as $n \to \infty$. For the order notation, you find the terms that will dominate as $m,\ n \to \infty$. For example, if the complexity is given by a function $F(n) = 3n^2 + 6n$, then we say $F(n) = \mathcal{O}(n^2)$. There are also many resources online for what the order notation means.*

## Extension Problem
With particularly large systems, often
the system of equations is instead solved using one of the following
three extremely common iterative methods:
- Jacobi method
- Gauss-Seidel method
- Successive Overrelaxation

Pick one of these methods and explain what the iterative scheme is and under what circumstances it converges. Give an example of a matrix for which the iterative scheme of your choice would not work. In what way is your method of choice better than row reduction?

*Commentary: This question should be done as a report. Pick one of the iterative methods and read about it and answer the questions I asked. You should submit a brief report which explains the method and summarizes when or when it doesn’t have convergence and if it’s better than row reduction.*

---

## Report
Let $Ax = b$\
Given an n x n real matrix $A$ (*coefficient matrix*) and a real n-vector $b$ (*right-hand side vector*), we need to find a vector $x$ (*unknown vector*).

$$
A_{n,n} = 
\begin{bmatrix}
a_{1,1} & a_{1,2} & \cdots & a_{1,n} \\
a_{2,1} & a_{2,2} & \cdots & a_{2,n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n,1} & a_{n,2} & \cdots & a_{n,n} 
\end{bmatrix}
$$

$A$ can be decomposed into a lower triangular part $L$, a diagonal part $D$, and an upper triangular part $U$.

$A = D + L + U$ where
$$
D = 
\begin{bmatrix}
a_{1,1} & 0 & \cdots & 0 \\
0 & a_{2,2} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & a_{n,n} 
\end{bmatrix}
{,\quad}
L = 
\begin{bmatrix}
0 & 0 & \cdots & 0 \\
a_{2,1} & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
a_{n,1} & a_{n,2} & \cdots & 0 
\end{bmatrix}
{,\quad}
U = 
\begin{bmatrix}
0 & a_{1,2} & \cdots & a_{1,n} \\
0 & 0 & \cdots & a_{2,n} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 0 
\end{bmatrix}
$$

### Iterative Scheme
#### Come up with an initial value for $x$
$x^{(0)}_1, x^{(0)}_2, \ldots, x^{(0)}_n$

#### Calculate the $k\text{-th}$ approximation / iteration of $x$
$$
\begin{aligned}
    \text{Component-wise Form: }\quad
    x^{(k+1)}_i &= \frac{1}{a_{ii}}\left( b_i - \sum_{j \neq i} a_{ij}x_{j}^{(k)} \right)
    , \quad i = 1, 2, \ldots , n

    \\

    \text{Vector Form: }\quad
    x^{(k+1)} &= D^{-1}\left(b - (L + U)x^{(k)}\right) \\
\end{aligned}
$$

#### Iterate until $\left\lVert b - Ax^{(k)} \right\rVert < \mathit{TOL}$ (Tolerance) or $\frac{\left\lVert x^{(k+1)} - x^{(k)} \right\rVert}{\left\lVert x^{(k+1)} \right\rVert}  < \mathit{TOL}$ (Tolerance)

### Convergence
Since $Ax = b$ and $A = D + L + U$\
From this we can rewrite the $k\text{-th}$ approximation / iteration of $x$ as
$$
\begin{aligned}
    (D + L + U)x &= b \\
    Dx + (L + U)x &= b \\
    Dx^{(k+1)} &= (L + U)x^{(k)} + b \\
    x^{(k+1)} &= -D^{-1}(L + U)x^{(k)} + D^{-1}b \\
    x^{(k+1)} &= Bx^{(k)} + c \\
\end{aligned}
$$
Where: $B = -D^{-1}(L + U)$ and $c = D^{-1}b$

For the Jacobi Method to converge, the spectral radius of the iteration matrix must be strictly less than 1. This means that the largest absolute value of the iteration matrix's eigenvalues must be strictly less than 1.

Let $B$ represent the iteration matrix and let $\lambda_1, \lambda_2, \ldots, \lambda_n$ be its corresponding eigenvalues.
$$
\begin{aligned}
    \rho \left(-D^{-1}(L + U) \right) &< 1 \\
    \rho \left(B \right) = \mathit{max} \left\{\left|\lambda_1 \right|, \ldots, \left|\lambda_n \right| \right\} &< 1
\end{aligned}
$$

The Jacobi Method guarentees a solution to the equation $Ax = b$ when $A$ is a strictly diagonally dominant matrix. This means that for the square matrix $A$ to be strictly diagonally dominant, for each row of the matrix, the magnitude of the diagonal entry must be larger than the sum of the magnitudes of the non-diagonal entries in that row
$$
\left| a_{ii} \right| > \sum_{j \neq i} \left| a_{ij} \right|
$$

---

### Additional Questions
An example of a matrix for which the iterative scheme of the Jacobi Method would not work is:

$$
A = 
\begin{bmatrix}
2 & 1 & 1 \\
1 & 2 & 3 \\
3 & 2 & 1 
\end{bmatrix}
$$

In addition to this matrix not being a strictly diagonally dominant matrix, the spectral radius of the iteration matrix $B$ is not strictly less than 1.

$$
B = -D^{-1}(L + U) = 
\begin{bmatrix}
-\frac{1}{2} & 0 & 0 \\
0 & -\frac{1}{2} & 0 \\
0 & 0 & -1 
\end{bmatrix}
\begin{bmatrix}
0 & 1 & 1 \\
1 & 0 & 3 \\
3 & 2 & 0 
\end{bmatrix}
=
\begin{bmatrix}
0 & -\frac{1}{2} & -\frac{1}{2} \\
-\frac{1}{2} & 0 & -\frac{3}{2} \\
-3 & -2 & 0 
\end{bmatrix}

\\

\rho \left(B \right) \approx 2.42563 > 1 \\
$$

Therefore the iteration matrix $B$ will not converge and the matrix $A$ would not work with the Jacobi Method.

---
#### The Jacobi Method is better than row reduction because
Since the $k\text{-th}$ approximation / iteration of $x$ can be written as

$$
\begin{aligned}
    \text{Component-wise Form: }\quad
    x^{(k+1)}_i &= \frac{1}{a_{ii}}\left( b_i - \sum_{j \neq i} a_{ij}x_{j}^{(k)} \right)
    , \quad i = 1, 2, \ldots , n

    \\

    \text{Vector Form: }\quad
    x^{(k+1)} &= D^{-1}\left(b - (L + U)x^{(k)}\right) \\
\end{aligned}
$$

The $k\text{-th}$ approximation / iteration of $x$ can be computed in parallel as the component-wise form can be split into parallel computations for $x^{(k+1)}_i$ as each $x^{(k+1)}_i$ does not rely on the previous value $x^{(k+1)}_{i - 1}$. This allows for this algorithm to be a good choice for large systems as its ability to be scaled with parallel computing as each component of the $k\text{-th}$ approximation / iteration of $x$ can be calculated in n-parallel threads instead of in a single thread. This also opens the possibility of applying this algorithm to GPUs which have hundreds or even thousands of smaller cores that perform several tasks at once and feature high memory bandwidth for large datasets.

In addition to this, the Jacobi Method is simpler than row reduction algorithm-wise as it calculates the $k\text{-th}$ approximation / iteration of $x$ until a certain condition is met or the algorithm reaches the iteration limit set compared to row reduction where partial pivoting may be needed to stabilize the system and it often has large floating point operations as a result of gaussian elimination.

$x^{(k+1)} = Bx^{(k)} + c$ \
Where: $B = -D^{-1}(L + U)$ and $c = D^{-1}b$

Iterate until $\left\lVert b - Ax^{(k)} \right\rVert < \mathit{TOL}$ (Tolerance), or \
$\frac{\left\lVert x^{(k+1)} - x^{(k)} \right\rVert}{\left\lVert x^{(k+1)} \right\rVert}  < \mathit{TOL}$ (Tolerance), or \
$k > \mathit{iteration\ limit}$



The Jacobi Method also is guaranteed to converge if the n x n and real matrix $A$ is strictly diagonally dominant from any initial guess $x$ for $x^{(0)}_1, x^{(0)}_2, \ldots, x^{(0)}_n$ and $Ax = b$

---

## References
[1] Yousef Saad. Iterative Methods for Sparse Linear Systems, Second Edition

[2] Ryan Reiff. Jacobi Iteration and Spectral Radius, A method for solving large linear systems of equations.

