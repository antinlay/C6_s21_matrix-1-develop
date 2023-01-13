#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int err_num = NO_OK;
  if (rows > 0 && columns > 0 && result) {
    result->rows = rows;
    result->columns = columns;
    result->matrix =
        (double **)calloc(rows, sizeof(double *) + columns * sizeof(double));
    if (result->matrix) {
      err_num = OK;
      for (int i = 0; i < rows; i++) {
        double *pntRows = (double *)(result->matrix + rows);
        result->matrix[i] = (pntRows + columns * i);
      }
    }
  }
  return err_num;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err_num = (!A || !B) ? NO_OK : OK;
  if (A->rows == B->rows && A->columns == B->columns && !err_num) {
    if (s21_create_matrix(A->rows, A->columns, result)) {
      err_num = NO_OK;
    } else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    }
  } else {
    err_num = ERR;
  }
  return err_num;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err_num = (!A || !B) ? NO_OK : OK;
  if (A->rows == B->rows && A->columns == B->columns && !err_num) {
    if (s21_create_matrix(A->rows, A->columns, result)) {
      err_num = NO_OK;
    } else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    }
  } else {
    err_num = ERR;
  }
  return err_num;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int code = 1;
  if (A->rows == B->rows && A->columns == B->columns && A && B) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= EPS) {
          code = 0;
          break;
        }
      }
      if (!code) break;
    }
  }
  return code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int err_num = !A ? NO_OK : OK;
  if (!err_num) {
    if (s21_create_matrix(A->rows, A->columns, result)) {
      err_num = NO_OK;
    } else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  }
  return err_num;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err_num = (!A || !B) ? NO_OK : OK;
  if (A->columns == B->rows && !err_num) {
    if (s21_create_matrix(A->rows, A->columns, result)) {
      err_num = NO_OK;
    } else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int k = 0; k < B->columns; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    }
  } else {
    err_num = ERR;
  }
  return err_num;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int err_num = !A ? NO_OK : OK;
  if (s21_create_matrix(A->columns, A->rows, result)) {
    err_num = NO_OK;
  } else {
    for (int i = 0; i < A->columns; i++) {
      for (int j = 0; j < A->rows; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  }
  return err_num;
}

int s21_determinant(matrix_t *A, double *result) {
  int err = 0, ok = 0;
  double quotient = 0, sign = 1;
  if (A == NULL || result == NULL) {
    err = NO_OK;
  } else if (A->rows == A->columns) {
    *result = 1;
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      matrix_t tmp = {0};
      s21_create_matrix(A->rows, A->columns, &tmp);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          tmp.matrix[i][j] = A->matrix[i][j];
        }
      }
      for (int i = 0; i < tmp.rows; i++) {
        if (tmp.matrix[i][i] == 0) {
          ok = swap_rows(&tmp, i);
          sign = -sign;
        }
        if (!ok) {
          for (int j = i + 1; j < tmp.rows; j++) {
            quotient = tmp.matrix[j][i] / tmp.matrix[i][i];
            for (int x = i; x < tmp.columns; x++) {
              tmp.matrix[j][x] = tmp.matrix[j][x] - quotient * tmp.matrix[i][x];
            }
          }
          *result *= tmp.matrix[i][i];
        } else {
          *result = 0;
          break;
        }
      }
      if (!ok) *result *= sign;
      s21_remove_matrix(&tmp);
    }
  } else {
    err = ERR;
  }
  return err;
}

void check_mtx(matrix_t *A, int x, int y, matrix_t *result) {
  int iDet, jDet;
  iDet = 0;
  for (int i = 0; i < (A->rows - 1); i++) {  // проверка индекса строки
    if (i == x) iDet = 1;
    jDet = 0;
    for (int j = 0; j < (A->rows - 1); j++) {  // проверка индекса столбца
      if (j == y) jDet = 1;
      result->matrix[i][j] = A->matrix[i + iDet][j + jDet];
    }
  }
}

int swap_rows(matrix_t *result, int index) {
  int code = OK;
  for (int i = index + 1; i < result->rows; i++) {
    if (!result->matrix[i][index]) {
      code = NO_OK;
    } else {
      code = OK;
      for (int j = 0; j < result->columns; j++) {
        double tmp;
        tmp = result->matrix[i][j];
        result->matrix[i][j] = result->matrix[index][j];
        result->matrix[index][j] = tmp;
      }
    }
    if (!code) break;
  }
  return code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int err = OK;
  if (A == NULL || result == NULL) {
    err = NO_OK;
  } else if (A->rows == A->columns) {
    matrix_t tmp = {0};
    double sign = 1;
    if (!s21_create_matrix(A->rows, A->columns, result)) {
      if (A->rows == 1) {
        result->matrix[0][0] = A->matrix[0][0];
      } else {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
            check_mtx(A, i, j, &tmp);
            s21_determinant(&tmp, &result->matrix[i][j]);
            s21_remove_matrix(&tmp);
            sign = pow(-1, i + j);
            result->matrix[i][j] = sign * result->matrix[i][j];
          }
        }
      }
    } else {
      err = NO_OK;
    }
  } else {
    err = ERR;
  }
  return err;
}

// int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
//   int err = OK;
//   double det = 0.0;
//   if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
//     err = NO_OK;
//   } else if (A->rows != A->columns) {
//     err = ERR;
//   } else {
//     err = s21_determinant(A, &det);
//     if (fabs(det) >= EPS && !err) {
//       matrix_t tmp = {0};
//       matrix_t buf = {0};
//       if (!err) {
//         err = s21_calc_complements(A, &tmp);
//       }
//       if (!err) {
//         err = s21_transpose(&tmp, &buf);
//       }
//       if (!err) {
//         err = s21_mult_number(&buf, (1.0 / det), result);
//       }
//       s21_remove_matrix(&buf);
//       s21_remove_matrix(&tmp);
//     } else {
//       err = ERR;
//     }
//   }
//   return err;
// }

double detr(double **mat, int dim) {
  double stat = 0;
  if (dim == 1) stat = mat[0][0];
  if (dim == 2) stat = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
  double det = 0.0, sign = 1.0;
  for (int i = 0; i < dim; i++) {
    matrix_t tmp = {0};
    s21_create_matrix(dim - 1, dim - 1, &tmp);
    for (int m = 1; m < dim; m++) {
      int p = 0;
      for (int n = 0; n < dim; n++) {
        if (n == i) continue;
        tmp.matrix[m - 1][p] = mat[m][n];
        p++;
      }
    }
    det += sign * mat[0][i] * detr(tmp.matrix, dim - 1);
    sign = -sign;
    s21_remove_matrix(&tmp);
  }
  if (stat == 0) stat = det;
  return stat;
}

void get_mtx(double **m, double **tmp, int skip_row, int skip_col, int size) {
  for (int row = 0, i = 0, j = 0; row < size; row++) {
    for (int col = 0; col < size; col++) {
      if (row != skip_row && col != skip_col) {
        tmp[i][j] = m[row][col];
        j++;
        if (j == size - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}

void inverse(matrix_t *A, matrix_t *result) {
  if (A->rows == 1) {
    result->matrix[0][0] = 1;
    return;
  }
  int size = A->rows;
  double **tmp = malloc(sizeof(double *) * size);
  for (int i = 0; i < size; i++) tmp[i] = malloc(sizeof(double) * size);
  int sign = 1;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      get_mtx(A->matrix, tmp, i, j, size);
      result->matrix[j][i] = sign * detr(tmp, size - 1);
      sign = -sign;
    }
  }
  for (int i = 0; i < size; i++) free(tmp[i]);
  free(tmp);
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int err = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    err = NO_OK;
  } else if (A->rows != A->columns) {
    err = ERR;
  } else {
    double det = 0.0;
    s21_determinant(A, &det);
    if (fabs(det) < EPS) err = ERR;
    matrix_t tmp = {0};
    s21_create_matrix(A->rows, A->columns, &tmp);
    inverse(A, &tmp);
    int size = A->rows;
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        result->matrix[i][j] = tmp.matrix[i][j] / det;
    s21_remove_matrix(&tmp);
  }
  return err;
}
