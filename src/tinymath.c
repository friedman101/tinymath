#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "tinymath.h"

int svd3(double A[3][3], double V[3][3], double s[3]) {
    gsl_matrix *A_gsl = gsl_matrix_alloc(3,3);
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++)
            gsl_matrix_set(A_gsl, i, j, A[i][j]);
    }

    gsl_matrix *V_gsl = gsl_matrix_alloc(3,3);
    gsl_vector *s_gsl = gsl_vector_alloc(3);
    gsl_vector *work = gsl_vector_alloc(3);
    
    int ret_val = gsl_linalg_SV_decomp(A_gsl, V_gsl, s_gsl, work);
    
    for (unsigned int i = 0; i < 3; i++) {
        s[i] = gsl_vector_get(s_gsl, i);
        for (unsigned int j = 0; j < 3; j++) {
            V[i][j] = gsl_matrix_get(V_gsl, i, j);
            A[i][j] = gsl_matrix_get(A_gsl, i, j);
        }
    }

    gsl_matrix_free(A_gsl);
    gsl_matrix_free(V_gsl);
    gsl_vector_free(s_gsl);
    gsl_vector_free(work);
    return ret_val;
}

void mat3_by_vec3(double y[3], double A[3][3], double x[3]) {
    for (unsigned int i = 0; i < 3; i++) {
        double sum = 0;
        for (unsigned int j = 0; j < 3; j++) {
            sum += A[i][j]*x[j];
        }
        y[i] = sum;
    }

}

void dcm2euler(double euler[3], double dcm[3][3]) {
    euler[0] = atan2(dcm[1][2], dcm[2][2]);
    euler[1] = asin(-dcm[0][2]);
    euler[2] = atan2(dcm[0][1], dcm[0][0]);
}

double uint8_std_mat(void *x_tmp, unsigned int row, unsigned int col, double mu) {
    uint8_t (*x)[col] = x_tmp;
    double var = 0;
    for (unsigned int i = 0; i < row; i++) {
        for (unsigned int j = 0; j < col; j++) {
            double tmp = ((double) x[i][j]) - mu;
            var += tmp*tmp;
        }
    }
    var = var/(row*col);
    return sqrt(var);
}

double uint8_mean_mat(void *x_tmp, unsigned int row, unsigned int col) {
    uint8_t (*x)[col] = x_tmp;
    double mu = 0;
    for (unsigned int i = 0; i < row; i++) {
        for (unsigned int j = 0; j < col; j++)
            mu += (double) x[i][j];
    }
    return mu/(row*col);
}

double dot3(double x1[3], double x2[3]) {
    double y = 0;
    for (unsigned int i = 0; i < 3; i++)
        y += x1[i]*x2[i];

    return y;
}

double mag3(double x[3]) {
    double y = dot3(x, x);
    return sqrt(y);
}

double ang_between(double x1[3], double x2[3]) {
    double mydot = dot3(x1, x2);
    double mag1 = mag3(x1);
    double mag2 = mag3(x2);

    return acos(mydot/mag1/mag2);
}


void normalize(double y[3], double x[3]) {
    double mymag = mag3(x);
    for (unsigned int i = 0; i < 3; i++)
        y[i] = y[i]/mymag;
}

double polyval(double x, double *p, unsigned int n) {
    double y = p[0];

    for (unsigned int i = 0; i < n; i++) {
        double x_pow = x;
        for (unsigned int j = 0; j < i; j++)
            x_pow = x_pow*x;
        y += p[i+1]*x_pow;
    }

    return y;
}

void cross(double y[3], double x1[3], double x2[3]) {
    y[0] = x1[1]*x2[2] - x1[2]*x2[1];
    y[1] = x1[2]*x2[0] - x1[0]*x2[2];
    y[2] = x1[0]*x2[1] - x1[1]*x2[0];
}

double triple_prod(double a[3], double b[3], double c[3]) {
    double tmp[3];
    cross(tmp, b, c);
    return dot3(a, tmp);
}

void radec2cart(double y[3], double RA, double dec) {
    y[0] = cos(dec)*cos(RA);
    y[1] = cos(dec)*sin(RA);
    y[2] = sin(dec);
}

void cart2radec(double *RA, double *dec, double x[3]) {
    *RA = atan2(x[1],x[0]);
    *dec = asin(x[2]);
}

int sign(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

void vec3_by_vec3_transpose(double Y[3][3], double u[3], double v[3]) {
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            Y[i][j] = u[i]*v[j];
        }
    }
}

void mat3_add(double Y[3][3], double U[3][3], double V[3][3]) {
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            Y[i][j] = U[i][j] + V[i][j];
        }
    }
}

void mat3_mult(double UV[3][3], double U[3][3], double V[3][3]) {
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            double sum = 0;
            for (unsigned int k = 0; k < 3; k++)
                sum += U[i][k]*V[k][j];
            UV[i][j] = sum;
        }
    }
}

void trans3(double Y[3][3]) {
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = i+1; j < 3; j++) {
            double tmp = Y[i][j];
            Y[i][j] = Y[j][i];
            Y[j][i] = tmp;
        }
    }
}

double det3(double X[3][3]) {
    double det = X[0][0]*((X[1][1]*X[2][2])
        - (X[2][1]*X[1][2])) -X[0][1]*(X[1][0]*X[2][2]
        - X[2][0]*X[1][2]) + X[0][2]*(X[1][0]*X[2][1]
        - X[2][0]*X[1][1]);
    return det;
}
