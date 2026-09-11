#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "utils.h"
#include "my_header.h"
static const double t4 = 0.45;


/*
 * Ridge regression predictor:
 *   y_hat = X (X^T X + λI)^(-1) X^T y
 *
 * Inputs:
 *   X      : pointer to double array (row-major) of size n * p
 *   y      : pointer to double array of length n
 *   n, p   : number of samples and predictors
 *   lambda : ridge penalty
 *
 * Output:
 *   y_hat  : pointer to double array of length n (predicted values)
 */
void ridge_predict_only(const double *X, const double *y,
                        int n, int p, double lambda,
                        double *y_hat)
{
    gsl_matrix_const_view Xv = gsl_matrix_const_view_array(X, n, p);
    gsl_vector_const_view yv = gsl_vector_const_view_array(y, n);
    gsl_vector_view yhatv = gsl_vector_view_array(y_hat, n);

    gsl_matrix *XtX = gsl_matrix_alloc(p, p);
    gsl_matrix *XtX_lambdaI = gsl_matrix_alloc(p, p);
    gsl_vector *Xty = gsl_vector_alloc(p);
    gsl_vector *beta = gsl_vector_alloc(p);

    // XtX = X^T * X
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Xv.matrix, &Xv.matrix, 0.0, XtX);

    // XtX + λI
    gsl_matrix_memcpy(XtX_lambdaI, XtX);
    for (int i = 0; i < p; i++) {
        double val = gsl_matrix_get(XtX_lambdaI, i, i) + lambda;
        gsl_matrix_set(XtX_lambdaI, i, i, val);
    }

    // Xty = X^T * y
    gsl_blas_dgemv(CblasTrans, 1.0, &Xv.matrix, &yv.vector, 0.0, Xty);

    // Solve (XtX + λI) * beta = Xty
    int signum;
    gsl_permutation *perm = gsl_permutation_alloc(p);
    gsl_linalg_LU_decomp(XtX_lambdaI, perm, &signum);
    gsl_linalg_LU_solve(XtX_lambdaI, perm, Xty, beta);

    // Compute y_hat = X * beta
    gsl_blas_dgemv(CblasNoTrans, 1.0, &Xv.matrix, beta, 0.0, &yhatv.vector);

    // Free memory
    gsl_matrix_free(XtX);
    gsl_matrix_free(XtX_lambdaI);
    gsl_vector_free(Xty);
    gsl_vector_free(beta);
    gsl_permutation_free(perm);
}
float generate_normal(const float sigma)
{

  // srand(1);
  float x, y, r2;

  do
  {
    /* choose x,y in uniform square (-1,-1) to (+1,+1) */
  //  x = -1 + 2 * ((double)rand() + 1.) / (1. + (double)RAND_MAX);
    //y = -1 + 2 * ((double)rand() + 1.) / (1. + (double)RAND_MAX);

// Clean, standard R-compatible uniform sampling between -1 and 1
x = -1.0 + 2.0 * unif_rand();
y = -1.0 + 2.0 * unif_rand();
    // printf("X=%2.5f \n",x);
    // printf("Y=%2.5f \n",y);
    /* see if it is in the unit circle */
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0);


  /* Box-Muller transform */
  return sigma * y * sqrt(-2.0 * log(r2) / r2);
}

/* The exponential distribution has the form

   p(x) dx = exp(-x/mu) dx/mu

   for x = 0 ... +infty */

double rexponential(const double mu)
{
  //double u = ((double)rand() + 1.) / (1. + (double)RAND_MAX);
  //return -mu * log1p(-u);
  return mu*exp_rand();
}

// Generae from truncated normal distribution

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(double a)
{
  // SAMPLER_DEBUG("ers_a_inf", a, R_PosInf);
  const double ainv = 1.0 / a;
  double x, z, rho;
  do
  {
    // x = rexp(ainv) + a; /* rexp works with 1/lambda */
    x = rexponential(ainv) + a;
    z= x - a;
    rho = exp(-0.5 * z * z);
  } while (unif_rand() > rho);
  return x;
}

/* Normal rejection sampling (a,inf) */
double nrs_a_inf(double a)
{
  // SAMPLER_DEBUG("nrs_a_inf", a, R_PosInf);
  // double x = -DBL_AX;
  double x = generate_normal(1.0);
  // double x = gsl_ran_ugaussian(r);
  while (x < a)
  {
    x = generate_normal(1.0);
    // x = gsl_ran_ugaussian(r);
  }
  return x;
}

double r_lefttruncnorm(double a, double mean, double sd)
{
  const double alpha = (a - mean) / sd;
  if (alpha < t4)
  {
    return mean + sd * nrs_a_inf(alpha);
  }
  else
  {
    return mean + sd * ers_a_inf(alpha);
  }
}
double r_righttruncnorm(double b, double mean, double sd)
{
  const double beta = (b - mean) / sd;
  /* Exploit symmetry: */
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0);
}





SEXP array_to_r_list(double **array, int rows, int *cols)
{
    // Allocate list of length = rows
    SEXP list = PROTECT(allocVector(VECSXP, rows));
    for (int i = 0; i < rows; i++)
    {
        // Allocate numeric vector of length = cols
        SEXP vec = PROTECT(allocVector(REALSXP, cols[i]));
        double *vec_ptr = REAL(vec);

        // Copy row i to numeric vector
        for (int j = 0; j < cols[i]; j++)
        {
            vec_ptr[j] = array[i][j];
        }
        // Set ith element of list
        SET_VECTOR_ELT(list, i, vec);
        UNPROTECT(1); // vec
    }
    UNPROTECT(1); // list
    return list;
}

// Function to convert a 2D C array to an R matrix
SEXP c_array_to_r_matrix_int(_Bool **array, int rows, int cols)
{
    // Step 1: Allocate a numeric vector of size rows * cols
    SEXP matrix = PROTECT(allocVector(INTSXP, rows * cols));
    int *matrix_data = INTEGER(matrix);

    // Step 2: Copy the data into the matrix in column-major order
    for (int col = 0; col < cols; ++col)
    {
        for (int row = 0; row < rows; ++row)
        {
            matrix_data[row + col * rows] = (int)array[row][col]; // R uses column-major order
        }
    }
    // Step 3: Set dimension attribute to convert it into a matrix
    SEXP dims = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dims)[0] = rows;
    INTEGER(dims)[1] = cols;
    setAttrib(matrix, R_DimSymbol, dims);
    UNPROTECT(2); // matrix, dims
    return matrix;
}

// Function to convert a 2D C array to an R matrix
SEXP c_array_to_r_matrix(double **array, int rows, int cols)
{
    // Step 1: Allocate a numeric vector of size rows * cols
    SEXP matrix = PROTECT(allocVector(REALSXP, rows * cols));
    double *matrix_data = REAL(matrix);

    // Step 2: Copy the data into the matrix in column-major order
    for (int col = 0; col < cols; ++col)
    {
        for (int row = 0; row < rows; ++row)
        {
            matrix_data[row + col * rows] = array[row][col]; // R uses column-major order
        }
    }

    // Step 3: Set dimension attribute to convert it into a matrix
    SEXP dims = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dims)
    [0] = rows;
    INTEGER(dims)
    [1] = cols;
    setAttrib(matrix, R_DimSymbol, dims);

    UNPROTECT(2); // matrix, dims
    return matrix;
}

double **r_list_vector_double_to_c(int listlength, SEXP LictVect)
{
    double **ListVect_c = malloc(listlength * sizeof(double *));
    for (int i = 0; i < listlength; i++)
    {
        SEXP mPM = VECTOR_ELT(LictVect, i);
        // Instead of getting dims via getAttrib, use the total length
        int sizeM = LENGTH(mPM);
        ListVect_c[i] = malloc(sizeM * sizeof(double));
        for (int j = 0; j < sizeM; j++)
        {
            ListVect_c[i][j] = REAL(mPM)[j];
        }
    }
    return ListVect_c;
}

double ***r_list_matrix_to_c(int listlength, SEXP ListMat)
{
    // ListMat is a list of matrices from R
    double ***newYY_arr = (double ***)malloc(listlength * sizeof(double **));

    for (int i = 0; i < listlength; i++)
    {
        SEXP mYY = VECTOR_ELT(ListMat, i);
        SEXP dimsYY = getAttrib(mYY, R_DimSymbol);
        //if (dimsYY == R_NilValue)
        //{
          //  UNPROTECT(2);
        //    error("ListMat element %d does not have dimension attributes.", i);
        //}
        int n_rows_yy = INTEGER(dimsYY)[0];
        int n_col_yy = INTEGER(dimsYY)[1];
        double *dataYY = REAL(mYY);
        // Allocate an array for row pointers if you intend to access newCC in 2D fashion.
        newYY_arr[i] = (double **)malloc(n_rows_yy * sizeof(double *));
        for (int r = 0; r < n_rows_yy; r++)
        {
            newYY_arr[i][r] = (double *)malloc(n_col_yy * sizeof(double));
            for (int c = 0; c < n_col_yy; c++)
            {
                newYY_arr[i][r][c] = dataYY[c * n_rows_yy + r]; // Still column-major indexing for newCC.
            }
        }
    }
    return newYY_arr;
}

_Bool ****r_list_list_matrix_to_c_bool(int listlength, SEXP ListListMat)
{
    // ListListMat is a list of list of matrices in R. Make sure that the data type is integer in R
    _Bool ****X0 = (_Bool ****)malloc(listlength * sizeof(_Bool ***));

    // Loop over each subgroup.
    for (int i = 0; i < listlength; i++)
    {
        SEXP subgroup = VECTOR_ELT(ListListMat, i);
        int n_platforms = LENGTH(subgroup);

        // Allocate space for the platforms within the subgroup.
        X0[i] = (_Bool ***)malloc(n_platforms * sizeof(_Bool **));

        // Loop over each platform.
        for (int j = 0; j < n_platforms; j++)
        {
            SEXP df = VECTOR_ELT(subgroup, j);

            // Get the matrix dimensions.
            SEXP dims = getAttrib(df, R_DimSymbol);
            int n_rows = INTEGER(dims)[0];
            int n_cols = INTEGER(dims)[1];

            // Get the pointer to the matrix data (numeric array).
            int *data_ptr = INTEGER(df);

            // Allocate an array to hold pointers to each row.
            X0[i][j] = (_Bool **)malloc(n_rows * sizeof(_Bool *));

            // Assign each row pointer; note R stores matrices in column-major order.
            for (int r = 0; r < n_rows; r++)
            {
                X0[i][j][r] = (_Bool *)malloc(n_cols * sizeof(_Bool));
                for (int c = 0; c < n_cols; c++)
                {
                    X0[i][j][r][c] = (_Bool)data_ptr[c * n_rows + r];
                }
            }
        }
    }
    return X0;
}

void free_r_list_list_matrix_to_c(double ****X0, int listlength, SEXP ListListMat)
{
    for (int i = 0; i < listlength; i++)
    {
        SEXP subgroup = VECTOR_ELT(ListListMat, i);
        int n_platforms = LENGTH(subgroup);

        for (int j = 0; j < n_platforms; j++)
        {
            SEXP df = VECTOR_ELT(subgroup, j);
            SEXP dims = getAttrib(df, R_DimSymbol);
            int n_rows = INTEGER(dims)[0];
            // printf("number of rows: %d \n", n_rows);

            for (int r = 0; r < n_rows; r++)
            {
                free(X0[i][j][r]);
                X0[i][j][r] = NULL;
            }
            free(X0[i][j]);
            X0[i][j] = NULL;
        }
        free(X0[i]);
        X0[i] = NULL;
    }
    free(X0);
    X0 = NULL;
}

double ****r_list_list_matrix_to_c(int listlength, SEXP ListListMat)
{
    // ListListMat is a list of list of matrices in R
    double ****X0 = (double ****)malloc(listlength * sizeof(double ***));

    // Loop over each subgroup.
    for (int i = 0; i < listlength; i++)
    {
        SEXP subgroup = VECTOR_ELT(ListListMat, i);
        int n_platforms = LENGTH(subgroup);

        // Allocate space for the platforms within the subgroup.
        X0[i] = (double ***)malloc(n_platforms * sizeof(double **));

        // Loop over each platform.
        for (int j = 0; j < n_platforms; j++)
        {
            SEXP df = VECTOR_ELT(subgroup, j);

            // If the matrix is not numeric but is logical, coerce it.
            if (!isReal(df) && isLogical(df))
            {
                df = coerceVector(df, REALSXP);
            }

            // Get the matrix dimensions.
            SEXP dims = getAttrib(df, R_DimSymbol);
            int n_rows = INTEGER(dims)[0];
            int n_cols = INTEGER(dims)[1];

            // Get the pointer to the matrix data (numeric array).
            double *data_ptr = REAL(df);

            // Allocate an array to hold pointers to each row.
            X0[i][j] = (double **)malloc(n_rows * sizeof(double *));

            // Assign each row pointer; note R stores matrices in column-major order.
            for (int r = 0; r < n_rows; r++)
            {
                X0[i][j][r] = (double *)malloc(n_cols * sizeof(double));
                for (int c = 0; c < n_cols; c++)
                {
                    X0[i][j][r][c] = data_ptr[c * n_rows + r];
                }
            }
        }
    }
    return X0;
}


static double mrf_log_weight(uint64_t state, int n_models,
                             double **theta, double nu)
{
    int selected = 0;
    double interaction = 0.0;
    for (int j = 0; j < n_models; ++j)
    {
        const unsigned int bit_j = (unsigned int)((state >> j) & UINT64_C(1));
        selected += (int)bit_j;
        interaction += bit_j * theta[j][j];
        for (int k = 0; k < j; ++k)
        {
            const unsigned int bit_k =
                (unsigned int)((state >> k) & UINT64_C(1));
            interaction += 2.0 * bit_j * bit_k * theta[j][k];
        }
    }
    return nu * selected + interaction;
}

void compute_mrf_log_normalizer(int n_models, double **theta, double nu,
                                double *log_normalizer)
{
    const uint64_t n_states = UINT64_C(1) << (unsigned int)n_models;
    double max_log_weight = -INFINITY;
    for (uint64_t state = 0; state < n_states; ++state)
    {
        max_log_weight = fmax(
            max_log_weight, mrf_log_weight(state, n_models, theta, nu));
    }
    double scaled_sum = 0.0;
    for (uint64_t state = 0; state < n_states; ++state)
    {
        scaled_sum += exp(
            mrf_log_weight(state, n_models, theta, nu) - max_log_weight);
    }
    *log_normalizer = max_log_weight + log(scaled_sum);
}

void sort_descending_index(int n, double *x, int *idx)
{
    int i, j;
    double a;
    int id;
    for (i = 0; i < n; i++)
        idx[i] = i;
    for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            if (x[i] <= x[j])
            {
                a = x[i];
                id = idx[i];
                idx[i] = idx[j];
                x[i] = x[j];
                idx[j] = id;
                x[j] = a;
            }
        }
    }
}


void mean_array_columns(int n, int n1, double **x, double *me)
{
    int i, l;
    for (i = 0; i < n1; i++)
    {
        me[i] = 0;
        for (l = 0; l < n; l++)
            me[i] += x[l][i] / n;
    }
}

_Bool bool_vectors_equal(int n, _Bool *u, _Bool *v)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (u[i] != v[i])
            return 0;
    }
    return 1;
}

double norm(int n, double *x)
{
    double normx = 0;
    int i;
    for (i = 0; i < n; i++)
    {
        normx += pow(x[i], 2);
        // printf("NormXXX==%f \n",x[i]);
    }
    return sqrt(normx);
}

double sum(int n, double *x)
{
    int i;
    double sum = 0;
    for (i = 0; i < n; i++)
        sum += x[i];
    return sum;
}

double mean(int n, double *x)
{
    int i;
    double me = 0;
    for (i = 0; i < n; i++)
        me += x[i];
    return me / n;
}
double var(int n, double *x)
{
    int i;
    double me = mean(n, x);
    double va = 0;
    for (i = 0; i < n; i++)
        va += (x[i] - me) * (x[i] - me);
    return va / (n - 1);
}

double *dvector(int nl, int nh)
{
    double *v;

    v = (double *)malloc((unsigned)(nh - nl + 1) * sizeof(double));
    if (!v)
        nrerror("allocation failure in dvector()");
    return v - nl;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;

    m = (double **)malloc((unsigned)(nrh - nrl + 1) * sizeof(double *));
    if (!m)
        nrerror("allocation failure 1 in dmatrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (double *)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
        if (!m[i])
            nrerror("allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }
    return m;
}

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    _Bool **m;

    m = (_Bool **)malloc((nrh - nrl + 1) * sizeof(_Bool *));
    if (!m)
        nrerror("allocation failure 1 in dmatrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++)
    {
        m[i] = (_Bool *)malloc((nch - ncl + 1) * sizeof(_Bool));
        if (!m[i])
            nrerror("allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }
    return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for (i = nrh; i >= nrl; i--)
        free((char *)(m[i] + ncl));

    free((char *)(m + nrl));
}

void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for (i = nrh; i >= nrl; i--)
        free((char *)(m[i] + ncl));

    free((char *)(m + nrl));
}

void nrerror(char error_text[])
{
    Rprintf("Utils run-time error...\n");
    Rprintf("%s\n", error_text);
    Rf_error("...now exiting to system...\n");
    // exit(1);
}
