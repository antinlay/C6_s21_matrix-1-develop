#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ERR 2
#define NO_OK 1
#define OK 0
#define EPS 1e-7

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);

int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

int s21_eq_matrix(matrix_t *A, matrix_t *B);

int swap_rows(matrix_t *result, int index);
void check_mtx(matrix_t *A, int x, int y, matrix_t *result);
double detr(double **mtx, int det);
void get_mtx(double **m, double **tmp, int skip_row, int skip_col, int size);
void inverse(matrix_t *A, matrix_t *result);

#endif  // SRC_S21_MATRIX_H_
