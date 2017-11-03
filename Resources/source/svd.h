// Functions to calculate the singular value decomposition of a matrix A
// and back-substitute to solve a matrix equation A x = b

// Singular value decomposition of matrix A of size m x n
int svdcmp(double **a, int m, int n, double *w, double **v);

// Back substitution of the singular value decomposition to solve for vector x
void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x);

// Pythagorean distance
double pythag(double a, double b);