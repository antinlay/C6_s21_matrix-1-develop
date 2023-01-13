#include <check.h>

#include "s21_matrix.h"

START_TEST(create_matrix) {
  matrix_t matrix1, matrix2;
  int ret = s21_create_matrix(0, 1, &matrix1);
  ck_assert_int_eq(ret, 1);
  ret = s21_create_matrix(3, 3, &matrix1);
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&matrix1);
  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&matrix1, &matrix2), 0);
  k = 0;
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&matrix1, &matrix2), 1);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);

  matrix_t wrong;
  ret = s21_create_matrix(-1, -2, &wrong);
  ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST(sum_sub_matrix) {
  matrix_t matrix1, matrix2, result;
  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  k = 0;
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  s21_sum_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_int_eq(result.matrix[i][j],
                       matrix2.matrix[i][j] + matrix1.matrix[i][j]);
    }
  }
  s21_remove_matrix(&result);
  s21_sub_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j],
                              matrix2.matrix[i][j] - matrix1.matrix[i][j],
                              1e-6);
    }
  }
  s21_remove_matrix(&result);

  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      srand(i + j);
      matrix1.matrix[i][j] = (double)rand() / (double)rand();
      srand((j + i) * 2);
      matrix2.matrix[i][j] = (double)rand() / (double)rand();
    }
  }

  s21_sum_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j],
                              matrix2.matrix[i][j] + matrix1.matrix[i][j],
                              1e-6);
    }
  }
  s21_remove_matrix(&result);

  s21_sub_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j],
                              matrix2.matrix[i][j] - matrix1.matrix[i][j],
                              1e-6);
    }
  }

  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);

  matrix_t wrong;
  s21_create_matrix(2, 3, &matrix1);
  s21_create_matrix(3, 2, &matrix2);
  int ret = s21_sub_matrix(&matrix1, &matrix2, &wrong);
  ck_assert_int_eq(ret, 2);
  ret = s21_sum_matrix(&matrix1, &matrix2, &wrong);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
}
END_TEST

START_TEST(mult_matrix) {
  matrix_t matrix1, matrix2, result;
  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  k = 0;
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  s21_mult_number(&matrix2, 2, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], matrix2.matrix[i][j] * 2,
                              1e-6);
    }
  }
  s21_remove_matrix(&result);
  int ret = s21_mult_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&result);

  matrix1.matrix[0][0] = 5.474748;
  matrix1.matrix[0][1] = -23.365346;
  matrix1.matrix[0][2] = 7.464645;
  matrix1.matrix[1][0] = 13.235363;
  matrix1.matrix[1][1] = -17.326774;
  matrix1.matrix[1][2] = -0.000034;
  matrix1.matrix[2][0] = -12.235567;
  matrix1.matrix[2][1] = 11.124526;
  matrix1.matrix[2][2] = 5.325634;

  matrix2.matrix[0][0] = 15.352534;
  matrix2.matrix[0][1] = 143.532636;
  matrix2.matrix[0][2] = 345.35356;
  matrix2.matrix[1][0] = 124.523552;
  matrix2.matrix[1][1] = -654.234562;
  matrix2.matrix[1][2] = 123.353578;
  matrix2.matrix[2][0] = -245.636465;
  matrix2.matrix[2][1] = 6324.235523;
  matrix2.matrix[2][2] = -2353.632542;
  s21_mult_matrix(&matrix1, &matrix2, &result);
  matrix_t calculated_mat;
  s21_create_matrix(3, 3, &calculated_mat);
  calculated_mat.matrix[0][0] = -4659.06324085756;
  calculated_mat.matrix[0][1] = 63280.3476665598;
  calculated_mat.matrix[0][2] = -18560.480561482;
  calculated_mat.matrix[1][0] = -1954.3776969092;
  calculated_mat.matrix[1][1] = 13235.232153486;
  calculated_mat.matrix[1][2] = 2433.63942021;
  calculated_mat.matrix[2][0] = -110.749147819;
  calculated_mat.matrix[2][1] = 24646.32859501;
  calculated_mat.matrix[2][2] = -15387.89622285;
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], calculated_mat.matrix[i][j],
                              1e-1);
    }
  }

  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);
  s21_remove_matrix(&calculated_mat);

  matrix_t wrong;
  s21_create_matrix(2, 3, &matrix1);
  s21_create_matrix(5, 6, &matrix2);
  ret = s21_mult_matrix(&matrix1, &matrix2, &wrong);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
}
END_TEST

START_TEST(transp_matrix) {
  matrix_t matrix1, matrix2, result;
  s21_create_matrix(3, 4, &matrix1);
  s21_create_matrix(4, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  matrix2.matrix[0][0] = 0;
  matrix2.matrix[0][1] = 4;
  matrix2.matrix[0][2] = 8;
  matrix2.matrix[1][0] = 1;
  matrix2.matrix[1][1] = 5;
  matrix2.matrix[1][2] = 9;
  matrix2.matrix[2][0] = 2;
  matrix2.matrix[2][1] = 6;
  matrix2.matrix[2][2] = 10;
  matrix2.matrix[3][0] = 3;
  matrix2.matrix[3][1] = 7;
  matrix2.matrix[3][2] = 11;
  s21_transpose(&matrix1, &result);
  for (int i = 0; i < result.rows; i++) {
    for (int j = 0; j < result.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], matrix2.matrix[i][j], 1e-6);
    }
  }
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(determ_matrix) {
  matrix_t matrix1;
  s21_create_matrix(3, 2, &matrix1);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  double res;
  int ret = s21_determinant(&matrix1, &res);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);

  s21_create_matrix(3, 3, &matrix1);
  k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  ret = s21_determinant(&matrix1, &res);
  ck_assert_int_eq(ret, 0);
  ck_assert_double_eq(res, 0);
  s21_remove_matrix(&matrix1);

  s21_create_matrix(3, 3, &matrix1);
  matrix1.matrix[0][0] = 0;
  matrix1.matrix[0][1] = 9;
  matrix1.matrix[0][2] = 5;

  matrix1.matrix[1][0] = 4;
  matrix1.matrix[1][1] = 3;
  matrix1.matrix[1][2] = -5;

  matrix1.matrix[2][0] = -1;
  matrix1.matrix[2][1] = 6;
  matrix1.matrix[2][2] = -4;
  ret = s21_determinant(&matrix1, &res);
  ck_assert_int_eq(ret, 0);
  ck_assert_double_eq(res, 324);
  s21_remove_matrix(&matrix1);

  double wrong;
  s21_create_matrix(2, 3, &matrix1);
  ret = s21_determinant(&matrix1, &wrong);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
}
END_TEST

START_TEST(calc_matrix) {
  matrix_t matrix1, result;
  s21_create_matrix(3, 2, &matrix1);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  int ret = s21_calc_complements(&matrix1, &result);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);

  s21_create_matrix(3, 3, &matrix1);
  k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  ret = s21_calc_complements(&matrix1, &result);
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&matrix1);

  s21_remove_matrix(&result);
  s21_create_matrix(3, 3, &matrix1);
  matrix1.matrix[0][0] = 1;
  matrix1.matrix[0][1] = 2;
  matrix1.matrix[0][2] = 3;

  matrix1.matrix[1][0] = 0;
  matrix1.matrix[1][1] = 4;
  matrix1.matrix[1][2] = 2;

  matrix1.matrix[2][0] = 5;
  matrix1.matrix[2][1] = 2;
  matrix1.matrix[2][2] = 1;
  s21_calc_complements(&matrix1, &result);
  matrix_t expected;
  s21_create_matrix(3, 3, &expected);
  expected.matrix[0][0] = 0;
  expected.matrix[0][1] = 10;
  expected.matrix[0][2] = -20;

  expected.matrix[1][0] = 4;
  expected.matrix[1][1] = -14;
  expected.matrix[1][2] = 8;

  expected.matrix[2][0] = -8;
  expected.matrix[2][1] = -2;
  expected.matrix[2][2] = 4;
  for (int i = 0; i < result.rows; i++) {
    for (int j = 0; j < result.columns; j++) {
      ck_assert_double_eq(result.matrix[i][j], expected.matrix[i][j]);
    }
  }
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&expected);

  matrix_t wrong;
  s21_create_matrix(2, 3, &matrix1);
  ret = s21_calc_complements(&matrix1, &wrong);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
}
END_TEST

START_TEST(inverse_matrix) {
  /* const int size = rand() % 100 + 1; */
  int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  m.matrix[0][0] = 2;
  m.matrix[0][1] = 5;
  m.matrix[0][2] = 7;
  m.matrix[1][0] = 6;
  m.matrix[1][1] = 3;
  m.matrix[1][2] = 4;
  m.matrix[2][0] = 5;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = -3;

  matrix_t res = {0};
  s21_create_matrix(size, size, &res);
  s21_inverse_matrix(&m, &res);

  matrix_t expected = {0};
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 1;
  expected.matrix[0][1] = -1;
  expected.matrix[0][2] = 1;
  expected.matrix[1][0] = -38;
  expected.matrix[1][1] = 41;
  expected.matrix[1][2] = -34;
  expected.matrix[2][0] = 27;
  expected.matrix[2][1] = -29;
  expected.matrix[2][2] = 24;

  ck_assert_int_eq(s21_eq_matrix(&expected, &res), 1);
  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&expected);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Core");
  TCase *matrix_test = tcase_create("Core");
  SRunner *sr = srunner_create(s1);
  suite_add_tcase(s1, matrix_test);

  tcase_add_test(matrix_test, create_matrix);
  tcase_add_test(matrix_test, sum_sub_matrix);
  tcase_add_test(matrix_test, mult_matrix);
  tcase_add_test(matrix_test, transp_matrix);
  tcase_add_test(matrix_test, determ_matrix);
  tcase_add_test(matrix_test, calc_matrix);
  tcase_add_test(matrix_test, inverse_matrix);

  srunner_run_all(sr, CK_VERBOSE);
  int err = srunner_ntests_failed(sr);
  srunner_free(sr);

  return err == 0 ? 0 : 1;
}
