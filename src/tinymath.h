#pragma once
#include <stdint.h>

int svd3(double A[3][3], double V[3][3], double s[3]);
void mat3_by_vec3(double y[3], double A[3][3], double x[3]);
void dcm2euler(double euler[3], double dcm[3][3]);
double uint8_std_mat(void *x_tmp, unsigned int row, unsigned int col, double mu);
double uint8_mean_mat(void *x_tmp, unsigned int row, unsigned int col);
double dot3(double x1[3], double x2[3]);
double mag3(double x[3]);
double ang_between(double x1[3], double x2[3]);
void normalize(double y[3], double x[3]);
double polyval(double x, double *p, unsigned int n);
void cross(double y[3], double x1[3], double x2[3]);
double triple_prod(double a[3], double b[3], double c[3]);
void radec2cart(double y[3], double RA, double dec);
void cart2radec(double *RA, double *dec, double x[3]);
int sign(double x);
void vec3_by_vec3_transpose(double Y[3][3], double u[3], double v[3]);
void mat3_add(double R[3][3], double U[3][3], double V[3][3]);
void mat3_mult(double UV[3][3], double U[3][3], double V[3][3]);
void trans3(double Y[3][3]);
double det3(double x[3][3]);
void vec3_by_scalar(double y[3], double x[3], double a);
