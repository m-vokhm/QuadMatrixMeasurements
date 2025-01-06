# A simple utility to estimate basic features of [QuadMatrix](https://github.com/m-vokhm/QuadMatrix).

Performs the most common operations on objects of the matrix classes from the QuadMatrix library, and for each combination of [type+operation+size] estimates the execution time (ms per operation) and accuracy (mean square error).

Examined classes:
- `com.mvohm.quadmatrix.BigDecimalMatrix`;
- `com.mvohm.quadmatrix.DoubleMatrix`;
- `com.mvohm.quadmatrix.QuadrupleMatrix`;
- `Jama.Matrix` *(https://math.nist.gov/javanumerics/jama/)*

Performed operations:
- Vector solution using LU-decomposition `(A * x = b)`
- Vector solution using LU-decomposition and iterative refinement (*except `JAMA`*) 
- Vector solution using Cholesky decomposition 
- Vector solution using Cholesky decomposition and iterative refinement (*except `JAMA`*)
- Matrix solution using LU-decomposition `(A * X = B)`
- Matrix solution using LU-decomposition and iterative refinement (*except `JAMA`*)
- Inversion
- Inversion with iterative refinement
- Multiplication

The current version of the code operates on matrices of sizes `50x50`, `100x100`, `200x200`, `400x400`.

The results are written to the file `./Results/stats_YYMMDD_HHMM.txt`, where `YYMMDD_HHMM` stands for local date and time, see an example in the Results directory.

Complete execution may take up to 48 hours or even more, depending on the performance of the processor. 




