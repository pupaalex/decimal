#include <check.h>
#include <stdlib.h>

#include "../src/s21_decimal.h"

///////////////////////////////////////////////////
// UPSACLE TESTS
///////////////////////////////////////////////////

START_TEST(UPSCALE_TEST1) {
  s21_decimal ds;
  ds.bits[0] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[2] = 0x00000007;
  ds.bits[3] = 0x00000000;
  for (int i = 0; i < 28; ++i) {
    int ret = s21_try_upscale(&ds);
    ck_assert(ret == 0);
  }
  int ret = s21_try_upscale(&ds);

  ck_assert(ret == 1);
}
END_TEST

START_TEST(UPSCALE_TEST2) {
  s21_decimal ds;
  ds.bits[0] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[2] = 0xffffffff;
  ds.bits[3] = 0x80000000;
  for (int i = 0; i < 28; ++i) {
    int ret = s21_try_upscale(&ds);
    ck_assert(ret == 1);
  }
  int ret = s21_try_upscale(&ds);

  ck_assert(ret == 1);
}
END_TEST

START_TEST(UPSCALE_TEST3) {
  s21_decimal ds;
  ds.bits[0] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[2] = 0x00000000;
  ds.bits[3] = 0x80000000;
  for (int i = 0; i < 28; ++i) {
    int ret = s21_try_upscale(&ds);
    ck_assert(ret == 0);
  }
  int ret = s21_try_upscale(&ds);

  ck_assert(ret == 1);

  ds.bits[0] = 0x1d14a021;
  ds.bits[1] = 0x9e548224;
  ds.bits[2] = 0x28000000;
  ds.bits[3] = 0x00000000;

  ret = s21_try_upscale(&ds);
  ck_assert(ret == 1);

  ds.bits[0] = 0x19999999;
  ds.bits[1] = 0x99999999;
  ds.bits[2] = 0x99999999;
  ds.bits[3] = 0x00000000;
  ret = s21_try_upscale(&ds);
  ck_assert(ret == 0);
  ds.bits[0] = 0x20000000;
  ds.bits[1] = 0x00000000;
  ds.bits[2] = 0x00000000;
  ds.bits[3] = 0x00000000;
  ret = s21_try_upscale(&ds);
  ck_assert(ret == 1);
}
END_TEST

Suite *dec_suite1(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal scale change tests");

  /* Core test case */
  tc_core = tcase_create("UPSCALE");

  tcase_add_test(tc_core, UPSCALE_TEST1);
  tcase_add_test(tc_core, UPSCALE_TEST2);
  tcase_add_test(tc_core, UPSCALE_TEST3);
  suite_add_tcase(s, tc_core);

  return s;
}
///////////////////////////////////////////////////
// ADD TESTS
///////////////////////////////////////////////////

START_TEST(ADD_TEST1) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000007;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = 0;
  ret = s21_add(ds, ds2, &res);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 9);
  ck_assert_int_eq(res.bits[3], 0);
  ck_assert_int_eq(ret, 0);
}
END_TEST

START_TEST(ADD_TEST2) {
  s21_decimal ds;
  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000001;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = 0;
  ret = s21_add(ds, ds2, &res);
  ck_assert_int_eq(ret, 1);
  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x80000000;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000001;
  ds2.bits[3] = 0x80000000;
  ret = s21_add(ds, ds2, &res);
  ck_assert_int_eq(ret, 2);
}
END_TEST

START_TEST(ADD_TEST3) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x0000000A;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x0000000A;
  ds2.bits[3] = 0x00010000;
  s21_decimal res1;
  s21_decimal res2;
  int ret = 0;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(res1.bits[2], 0);
  ck_assert_int_eq(res1.bits[1], 0);
  ck_assert_int_eq(res1.bits[0], 11);
  ck_assert_int_eq(res1.bits[3], 0);
  ck_assert_int_eq(ret, 0);

  ret = s21_add(ds2, ds, &res2);
  ck_assert_int_eq(res2.bits[2], res1.bits[2]);
  ck_assert_int_eq(res2.bits[1], res1.bits[1]);
  ck_assert_int_eq(res2.bits[0], res1.bits[0]);
  ck_assert_int_eq(res2.bits[3], res1.bits[3]);
  ck_assert_int_eq(ret, 0);
}
END_TEST

START_TEST(ADD_TEST4) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000009;
  ds.bits[3] = 0x00010000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000001;
  ds2.bits[3] = 0x00000000;
  s21_decimal res1, res2;
  int ret = 0;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(res1.bits[2], 0);
  ck_assert_int_eq(res1.bits[1], 0);
  ck_assert_int_eq(res1.bits[0], 19);
  ck_assert_int_eq(res1.bits[3], 0x00010000);

  ret = s21_add(ds2, ds, &res2);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(res2.bits[2], 0);
  ck_assert_int_eq(res2.bits[1], 0);
  ck_assert_int_eq(res2.bits[0], 19);
  ck_assert_int_eq(res2.bits[3], 0x00010000);
}
END_TEST

START_TEST(ADD_TEST5) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000308;
  ds.bits[0] = 0x1a263555;
  ds.bits[3] = 0x000b0000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000001;
  ds2.bits[3] = 0x000c0000;
  s21_decimal res;
  int ret = 0;
  ret = s21_add(ds, ds2, &res);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0x00001e51);
  ck_assert_int_eq(res.bits[0], 0x057e1553);
  ck_assert_int_eq(res.bits[3], 0x000c0000);
  ck_assert_int_eq(ret, 0);

  ds.bits[2] = 0x6bb4afe4;
  ds.bits[1] = 0xcf2607ee;
  ds.bits[0] = 0x35555555;
  ds.bits[3] = 0x001b0000;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000001;
  ds2.bits[3] = 0x001c0000;
  ret = s21_add(ds, ds2, &res);
  ck_assert_int_eq(res.bits[3], 0x001b0000);
  ck_assert_int_eq(res.bits[2], 0x6bb4afe4);
  ck_assert_int_eq(res.bits[1], 0xcf2607ee);
  ck_assert_int_eq(res.bits[0], 0x35555555);
  ck_assert_int_eq(ret, 0);
}
END_TEST

START_TEST(ADD_TEST6) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x0000000A;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x0000000A;
  ds2.bits[3] = 0x80000000;
  s21_decimal res1;
  int ret = 0;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(res1.bits[2], 0);
  ck_assert_int_eq(res1.bits[1], 0);
  ck_assert_int_eq(res1.bits[0], 0);
  ck_assert_int_eq(res1.bits[3], 0);
  ck_assert_int_eq(ret, 0);

  ds2.bits[2] = 0xffffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x80000000;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(res1.bits[2], 0xffffffff);
  ck_assert_int_eq(res1.bits[1], 0xffffffff);
  ck_assert_int_eq(res1.bits[0], 0xfffffff5);
  ck_assert_int_eq(res1.bits[3], 0x80000000);
  ck_assert_int_eq(ret, 0);
}
END_TEST

START_TEST(ADD_TEST7) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x001c0000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x801c0000;
  s21_decimal res1;
  s21_set_zeros(&res1);
  int ret = 0;
  for (int i = 0; i < 100; ++i) {
    ret = s21_add(res1, ds, &res1);
    ck_assert(ret == 0);
  }
  // printf("%x\n%x\n%x\n%x\n", res1.bits[2], res1.bits[1], res1.bits[0],
  // res1.bits[3]);
  for (int i = 0; i < 99; ++i) {
    ret = s21_add(res1, ds2, &res1);
    ck_assert(ret == 0);
  }
  // printf("%x\n%x\n%x\n%x\n", res1.bits[2], res1.bits[1], res1.bits[0],
  // res1.bits[3]);
  ck_assert_int_eq(res1.bits[2], 0x00000000);
  ck_assert_int_eq(res1.bits[1], 0xffffffff);
  ck_assert_int_eq(res1.bits[0], 0xffffffff);
  ck_assert_int_eq(res1.bits[3], 0x001c0000);
  ret = s21_add(res1, ds2, &res1);
  ck_assert(ret == 0);
  ck_assert_int_eq(res1.bits[2], 0x00000000);
  ck_assert_int_eq(res1.bits[1], 0x00000000);
  ck_assert_int_eq(res1.bits[0], 0x00000000);
  ck_assert_int_eq(res1.bits[3], 0x00000000);
  ret = s21_add(res1, ds2, &res1);
  ck_assert(ret == 0);
  ck_assert_int_eq(res1.bits[2], 0x00000000);
  ck_assert_int_eq(res1.bits[1], 0xffffffff);
  ck_assert_int_eq(res1.bits[0], 0xffffffff);
  ck_assert_int_eq(res1.bits[3], 0x801c0000);
  ck_assert_int_eq(ret, 0);
}
END_TEST

START_TEST(ADD_TEST8) {
  s21_decimal ds;
  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x001c0000;
  s21_decimal ds2;
  ds2.bits[2] = 0xffffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x001c0000;
  s21_decimal res1;
  int ret = 0;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(res1.bits[2], 0x33333333);
  ck_assert_int_eq(res1.bits[1], 0x33333333);
  ck_assert_int_eq(res1.bits[0], 0x33333333);
  ck_assert_int_eq(res1.bits[3], 0x001b0000);

  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x801c0000;

  ds2.bits[2] = 0xffffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x801c0000;
  ret = 0;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(res1.bits[2], 0x33333333);
  ck_assert_int_eq(res1.bits[1], 0x33333333);
  ck_assert_int_eq(res1.bits[0], 0x33333333);
  ck_assert_int_eq(res1.bits[3], 0x801b0000);
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000005;
  ds.bits[3] = 0x00100000;

  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x0000004b;
  ds2.bits[3] = 0x80110000;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(res1.bits[2], 0x00000000);
  ck_assert_int_eq(res1.bits[1], 0x00000000);
  ck_assert_int_eq(res1.bits[0], 0x00000019);
  ck_assert_int_eq(res1.bits[3], 0x80110000);
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000005;
  ds.bits[3] = 0x80100000;

  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x0000004b;
  ds2.bits[3] = 0x00110000;
  ret = s21_add(ds, ds2, &res1);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(res1.bits[2], 0x00000000);
  ck_assert_int_eq(res1.bits[1], 0x00000000);
  ck_assert_int_eq(res1.bits[0], 0x00000019);
  ck_assert_int_eq(res1.bits[3], 0x00110000);

  ds.bits[2] = 0x204fce57;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x001c0000;

  ds2.bits[2] = 0xefffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x00000000;
  ret = s21_add(ds, ds2, &res1);
  ck_assert(ret == 0);
  ck_assert_int_eq(res1.bits[2], 0xf0000000);
  ck_assert_int_eq(res1.bits[1], 0x00000000);
  ck_assert_int_eq(res1.bits[0], 0x00000000);
  ck_assert_int_eq(res1.bits[3], 0x00000000);

  ds.bits[2] = 0x204fce57;
  ds.bits[1] = 0x00000005;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x001c0000;

  ds2.bits[2] = 0xefffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x00000000;
  ret = s21_add(ds, ds2, &res1);

  ck_assert_int_eq(res1.bits[2], 0xf0000000);
  ck_assert_int_eq(res1.bits[1], 0x00000000);
  ck_assert_int_eq(res1.bits[0], 0x00000000);
  ck_assert_int_eq(res1.bits[3], 0x00000000);

  ck_assert_int_eq(ret, 0);
}
END_TEST

Suite *dec_suite2(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal add test");
  tc_core = tcase_create("ADD");
  tcase_add_test(tc_core, ADD_TEST1);
  tcase_add_test(tc_core, ADD_TEST2);
  tcase_add_test(tc_core, ADD_TEST3);
  tcase_add_test(tc_core, ADD_TEST4);
  tcase_add_test(tc_core, ADD_TEST5);
  tcase_add_test(tc_core, ADD_TEST6);
  tcase_add_test(tc_core, ADD_TEST7);
  tcase_add_test(tc_core, ADD_TEST8);
  suite_add_tcase(s, tc_core);

  return s;
}

///////////////////////////////////////////////////
// EQUAL TESTS
///////////////////////////////////////////////////

START_TEST(EQUAL_TEST1) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000000;
  ds.bits[3] = 0x00010000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x00000000;
  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));
}
END_TEST

START_TEST(EQUAL_TEST2) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000002;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000004;
  ds2.bits[3] = 0x00000000;
  ck_assert(s21_is_less_or_equal(ds, ds2));
}
END_TEST

START_TEST(EQUAL_TEST3) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000000;
  ds.bits[3] = 0x800a0000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x001c0000;

  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));

  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x80000000;

  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));

  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000000;
  ds.bits[3] = 0x000a0000;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x001c0000;

  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));

  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000000;
  ds.bits[3] = 0x000a0000;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x801c0000;

  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));
}
END_TEST

START_TEST(EQUAL_TEST4) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000001;
  ds.bits[0] = 0xffffffa4;
  ds.bits[3] = 0x00020000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x3333332a;
  ds2.bits[3] = 0x00010000;
  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));

  ds.bits[3] = 0x80020000;

  ck_assert(s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(!s21_is_greater_or_equal(ds, ds2));
  ck_assert(!s21_is_equal(ds, ds2));
  ck_assert(s21_is_not_equal(ds, ds2));
  ds2.bits[3] = 0x80010000;

  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));
  ds.bits[3] = 0x00020000;

  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(!s21_is_less_or_equal(ds, ds2));
  ck_assert(s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(!s21_is_equal(ds, ds2));
  ck_assert(s21_is_not_equal(ds, ds2));
}
END_TEST

START_TEST(EQUAL_TEST5) {
  s21_decimal ds;
  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0xffffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x00000000;
  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(s21_is_equal(ds, ds2));
  ck_assert(!s21_is_not_equal(ds, ds2));
  ds2.bits[3] = 0x80000000;
  ck_assert(!s21_is_less(ds, ds2));
  ck_assert(!s21_is_less_or_equal(ds, ds2));
  ck_assert(s21_is_greater(ds, ds2));
  ck_assert(s21_is_greater_or_equal(ds, ds2));
  ck_assert(!s21_is_equal(ds, ds2));
  ck_assert(s21_is_not_equal(ds, ds2));

  ds.bits[2] = 0x7fffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x001c0000;
  ds2.bits[2] = 0xffffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xfffffff1;
  ds2.bits[3] = 0x001c0000;
  ck_assert(s21_is_less(ds, ds2));
  ck_assert(s21_is_less_or_equal(ds, ds2));
  ck_assert(!s21_is_greater(ds, ds2));
  ck_assert(!s21_is_greater_or_equal(ds, ds2));
  ck_assert(!s21_is_equal(ds, ds2));
  ck_assert(s21_is_not_equal(ds, ds2));
}
END_TEST

Suite *dec_suite3(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal equal test");
  tc_core = tcase_create("EQUAL");
  tcase_add_test(tc_core, EQUAL_TEST1);
  tcase_add_test(tc_core, EQUAL_TEST2);
  tcase_add_test(tc_core, EQUAL_TEST3);
  tcase_add_test(tc_core, EQUAL_TEST4);
  tcase_add_test(tc_core, EQUAL_TEST5);
  suite_add_tcase(s, tc_core);

  return s;
}

///////////////////////////////////////////////////
// SUB TESTS
///////////////////////////////////////////////////

START_TEST(SUB_TEST1) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000000;
  ds.bits[3] = 0x00010000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000002;
  ds.bits[3] = 0x001b0000;

  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x0000000a;
  ds2.bits[3] = 0x001c0000;

  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 1);
  ck_assert_int_eq(res.bits[3], 0x001b0000);
}
END_TEST

START_TEST(SUB_TEST2) {
  s21_decimal ds;
  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0xffffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = 0;
  ret = s21_sub(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ds.bits[3] = 0x80000000;
  ds2.bits[3] = 0x80000000;

  ret = s21_sub(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0x80000000);

  ds2.bits[3] = 0x00000000;
  ret = s21_sub(ds, ds2, &res);
  ck_assert(ret == 2);
  ds2.bits[3] = 0x80000000;
  ds.bits[3] = 0x00000000;
  ret = s21_sub(ds, ds2, &res);
  ck_assert(ret == 1);
}
END_TEST

START_TEST(SUB_TEST3) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000000;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x80000000;
  s21_decimal res;
  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x000000c7;
  ds.bits[0] = 0xffffe962;
  ds.bits[3] = 0x00020000;

  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x0000306a;
  ds2.bits[0] = 0xd00d884c;
  ds2.bits[3] = 0x00080000;

  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0x0beb9193);
  ck_assert_int_eq(res.bits[0], 0xd6d61434);
  ck_assert_int_eq(res.bits[3], 0x00080000);
  ds.bits[3] = 0x80020000;
  ds2.bits[3] = 0x80080000;
  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0x0beb9193);
  ck_assert_int_eq(res.bits[0], 0xd6d61434);
  ck_assert_int_eq(res.bits[3], 0x80080000);

  ds.bits[3] = 0x00020000;

  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0x0bebf269);
  ck_assert_int_eq(res.bits[0], 0x76f124cc);
  ck_assert_int_eq(res.bits[3], 0x00080000);

  ds.bits[3] = 0x80020000;
  ds2.bits[3] = 0x00080000;

  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0x0bebf269);
  ck_assert_int_eq(res.bits[0], 0x76f124cc);
  ck_assert_int_eq(res.bits[3], 0x80080000);

  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x00000000;

  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000006;
  ds2.bits[3] = 0x00010000;
  ck_assert(s21_sub(ds, ds2, &res) == 0);
  ck_assert_int_eq(res.bits[2], 0xffffffff);
  ck_assert_int_eq(res.bits[1], 0xffffffff);
  ck_assert_int_eq(res.bits[0], 0xfffffffe);
  ck_assert_int_eq(res.bits[3], 0x00000000);
}
END_TEST

Suite *dec_suite4(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal SUB test");
  tc_core = tcase_create("SUB");
  tcase_add_test(tc_core, SUB_TEST1);
  tcase_add_test(tc_core, SUB_TEST2);
  tcase_add_test(tc_core, SUB_TEST3);
  suite_add_tcase(s, tc_core);

  return s;
}

///////////////////////////////////////////////////
// MUL TESTS
///////////////////////////////////////////////////

START_TEST(MUL_TEST1) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000001;
  ds.bits[3] = 0x00010000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = s21_mul(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 2);
  ck_assert_int_eq(res.bits[3], 0x00010000);
}
END_TEST

START_TEST(MUL_TEST2) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x80000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = s21_mul(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0x00000001);
  ck_assert_int_eq(res.bits[0], 0xfffffffe);
  ck_assert_int_eq(res.bits[3], 0x80000000);
  ds.bits[2] = 0xffffffff;
  ret = s21_mul(ds, ds2, &res);
  ck_assert(ret == 2);
}
END_TEST

START_TEST(MUL_TEST3) {
  s21_decimal ds;
  ds.bits[2] = 0x0bffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x13ffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = s21_mul(ds, ds2, &res);
  ck_assert(ret == 1);
  ds.bits[3] = 0x001c0000;
  ds2.bits[3] = 0x001c0000;
  ret = s21_mul(ds, ds2, &res);
  ck_assert(ret == 0);
}
END_TEST

Suite *dec_suite5(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal MUL test");
  tc_core = tcase_create("MULTIP");
  tcase_add_test(tc_core, MUL_TEST1);
  tcase_add_test(tc_core, MUL_TEST2);
  tcase_add_test(tc_core, MUL_TEST3);
  // tcase_add_test(tc_core, SUB_TEST2);
  // tcase_add_test(tc_core, SUB_TEST3);
  // tcase_add_test(tc_core, SUB_TEST4);
  suite_add_tcase(s, tc_core);

  return s;
}

///////////////////////////////////////////////////
// DIV TESTS
///////////////////////////////////////////////////

START_TEST(DIV_TEST1) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000002;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 1);
  ck_assert_int_eq(res.bits[3], 0x00000000);
}
END_TEST

START_TEST(DIV_TEST2) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000004;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0x00000002);
  ck_assert_int_eq(res.bits[3], 0x00000000);
}
END_TEST

START_TEST(DIV_TEST3) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000003;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x0000000f;
  ds2.bits[3] = 0x00010000;

  s21_decimal res;
  int ret = 0;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0x00000002);
  ck_assert_int_eq(res.bits[3], 0x00000000);

  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000040;
  ds.bits[3] = 0x00000000;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0x00000020);
  ck_assert_int_eq(res.bits[3], 0x00000000);

  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000003;
  ds.bits[3] = 0x00000000;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0x0000000f);
  ck_assert_int_eq(res.bits[3], 0x00010000);
}
END_TEST

START_TEST(DIV_TEST4) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000000;
  ds.bits[0] = 0x00000064;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000003;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0x6bb4afe4);
  ck_assert_int_eq(res.bits[1], 0xcf2607ee);
  ck_assert_int_eq(res.bits[0], 0x35555555);
  ck_assert_int_eq(res.bits[3], 0x001b0000);
}
END_TEST

START_TEST(DIV_TEST5) {
  s21_decimal ds;
  ds.bits[2] = 0xd7695fc9;
  ds.bits[1] = 0x9e4c0fdc;
  ds.bits[0] = 0x6aaaaaab;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000002;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0x6bb4afe4);
  ck_assert_int_eq(res.bits[1], 0xcf2607ee);
  ck_assert_int_eq(res.bits[0], 0x35555556);
  ck_assert_int_eq(res.bits[3], 0x00000000);
}
END_TEST

START_TEST(DIV_TEST6) {
  s21_decimal ds;
  ds.bits[2] = 0x00000000;
  ds.bits[1] = 0x00000005;
  ds.bits[0] = 0x00000000;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000003;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = 0;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0xe74be507);
  ck_assert_int_eq(res.bits[1], 0x9082aaaa);
  ck_assert_int_eq(res.bits[0], 0xaaaaaaab);
  ck_assert_int_eq(res.bits[3], 0x00130000);

  ds.bits[3] = 0x80010000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0xe74be507);
  ck_assert_int_eq(res.bits[1], 0x9082aaaa);
  ck_assert_int_eq(res.bits[0], 0xaaaaaaab);
  ck_assert_int_eq(res.bits[3], 0x80140000);

  ds2.bits[3] = 0x80000000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0xe74be507);
  ck_assert_int_eq(res.bits[1], 0x9082aaaa);
  ck_assert_int_eq(res.bits[0], 0xaaaaaaab);
  ck_assert_int_eq(res.bits[3], 0x00140000);

  ds.bits[3] = 0x00010000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0xe74be507);
  ck_assert_int_eq(res.bits[1], 0x9082aaaa);
  ck_assert_int_eq(res.bits[0], 0xaaaaaaab);
  ck_assert_int_eq(res.bits[3], 0x80140000);

  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x001c0000;
  ds2.bits[2] = 0xffffffff;
  ds2.bits[1] = 0xffffffff;
  ds2.bits[0] = 0xffffffff;
  ds2.bits[3] = 0x00000000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0x00000000);
  ck_assert_int_eq(res.bits[1], 0x00000000);
  ck_assert_int_eq(res.bits[0], 0x00000001);
  ck_assert_int_eq(res.bits[3], 0x001c0000);

  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xfffffff1;
  ds.bits[3] = 0x001c0000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0x00000000);
  ck_assert_int_eq(res.bits[1], 0x00000000);
  ck_assert_int_eq(res.bits[0], 0x00000001);
  ck_assert_int_eq(res.bits[3], 0x001c0000);

  ds.bits[3] = 0x801c0000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0x00000000);
  ck_assert_int_eq(res.bits[1], 0x00000000);
  ck_assert_int_eq(res.bits[0], 0x00000001);
  ck_assert_int_eq(res.bits[3], 0x801c0000);

  ds2.bits[3] = 0x80000000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0x00000000);
  ck_assert_int_eq(res.bits[1], 0x00000000);
  ck_assert_int_eq(res.bits[0], 0x00000001);
  ck_assert_int_eq(res.bits[3], 0x001c0000);

  ds.bits[3] = 0x001c0000;

  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0x00000000);
  ck_assert_int_eq(res.bits[1], 0x00000000);
  ck_assert_int_eq(res.bits[0], 0x00000001);
  ck_assert_int_eq(res.bits[3], 0x801c0000);
}
END_TEST

START_TEST(DIV_TEST7) {
  s21_decimal ds;
  ds.bits[2] = 0x0;
  ds.bits[1] = 0x0;
  ds.bits[0] = 0x0;
  ds.bits[3] = 0x00000000;
  s21_decimal ds2;
  ds2.bits[2] = 0x00000000;
  ds2.bits[1] = 0x00000000;
  ds2.bits[0] = 0x00000000;
  ds2.bits[3] = 0x00000000;
  s21_decimal res;
  int ret = 0;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 3);
  ds.bits[3] = 0x001c0000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 3);
  ds2.bits[3] = 0x001c0000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 3);
  ds2.bits[3] = 0x801c0000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 3);
  ds.bits[3] = 0x801c0000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 3);
  ds2.bits[0] = 0x0000001;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);
  ds.bits[3] = 0x001c0000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0x80000000);
  ds2.bits[3] = 0x001c0000;
  ret = s21_div(ds, ds2, &res);
  ck_assert(ret == 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0x00000000);
}
END_TEST

Suite *dec_suite6(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal DIV test");
  tc_core = tcase_create("DIVISION");
  tcase_add_test(tc_core, DIV_TEST1);
  tcase_add_test(tc_core, DIV_TEST2);
  tcase_add_test(tc_core, DIV_TEST3);
  tcase_add_test(tc_core, DIV_TEST4);
  tcase_add_test(tc_core, DIV_TEST5);
  tcase_add_test(tc_core, DIV_TEST6);
  tcase_add_test(tc_core, DIV_TEST7);
  suite_add_tcase(s, tc_core);

  return s;
}

///
// FROM NUM TO DEC
///

START_TEST(CONV_TEST1) {
  s21_decimal ds;
  int ret = s21_from_float_to_decimal(1.0, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000001);
  ck_assert_int_eq(ds.bits[3], 0x00000000);
  ret = s21_from_float_to_decimal(-1.0, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000001);
  ck_assert_int_eq(ds.bits[3], 0x80000000);
  float small = 0.0000000000000000000000000001;
  ret = s21_from_float_to_decimal(small, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000001);
  ck_assert_int_eq(ds.bits[3], 0x001c0000);
  ret = s21_from_float_to_decimal(-small, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000001);
  ck_assert_int_eq(ds.bits[3], 0x801c0000);
}
END_TEST

START_TEST(CONV_TEST2) {
  s21_decimal ds;
  float extrasmall = 0.000000000000000000000000000099;
  int ret = 0;
  ret = s21_from_float_to_decimal(extrasmall, &ds);
  ck_assert(ret == 1);
  ret = s21_from_float_to_decimal(extrasmall, &ds);
  ck_assert(ret == 1);
  float big = 7922816251426433759354395033.0;
  ret = s21_from_float_to_decimal(big, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x0098967f);
  ck_assert_int_eq(ds.bits[3], 0x00000000);
  ret = s21_from_float_to_decimal(-big, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x0098967f);
  ck_assert_int_eq(ds.bits[3], 0x80000000);
  float extrabig = 79228162514264337593543950336.0;
  ret = s21_from_float_to_decimal(extrabig, &ds);
  ck_assert(ret == 1);
  ret = s21_from_float_to_decimal(-extrabig, &ds);
  ck_assert(ret == 1);
}
END_TEST

START_TEST(CONV_TEST3) {
  s21_decimal ds;
  float a = 5.5;
  int ret = 0;
  ret = s21_from_float_to_decimal(a, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000037);
  ck_assert_int_eq(ds.bits[3], 0x00010000);
  ret = s21_from_float_to_decimal(-a, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000037);
  ck_assert_int_eq(ds.bits[3], 0x80010000);
  a = 0.5;
  ret = s21_from_float_to_decimal(a, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000005);
  ck_assert_int_eq(ds.bits[3], 0x00010000);
  ret = s21_from_float_to_decimal(-a, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00000005);
  ck_assert_int_eq(ds.bits[3], 0x80010000);
  a = 545.342432432;
  ret = s21_from_float_to_decimal(a, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00533670);
  ck_assert_int_eq(ds.bits[3], 0x00040000);
  ret = s21_from_float_to_decimal(-a, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x00533670);
  ck_assert_int_eq(ds.bits[3], 0x80040000);

  ck_assert(s21_set_scale(30, &ds));

  int t = 0x7fffffff;
  ret = s21_from_int_to_decimal(t, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x7fffffff);
  ck_assert_int_eq(ds.bits[3], 0x00000000);
  ret = s21_from_int_to_decimal(-t, &ds);
  ck_assert(ret == 0);
  ck_assert_int_eq(ds.bits[2], 0x00000000);
  ck_assert_int_eq(ds.bits[1], 0x00000000);
  ck_assert_int_eq(ds.bits[0], 0x7fffffff);
  ck_assert_int_eq(ds.bits[3], 0x80000000);
}
END_TEST

START_TEST(CONV_TEST4) {
  s21_decimal ds = {0};
  int ans = 0;
  ck_assert(s21_from_decimal_to_int(ds, &ans) == 0);
  ck_assert_int_eq(ans, 0);
  ds.bits[3] = 0x80000000;
  ck_assert(s21_from_decimal_to_int(ds, &ans) == 0);
  ck_assert_int_eq(ans, 0);
  ds.bits[3] = 0x801c0000;
  ck_assert(s21_from_decimal_to_int(ds, &ans) == 0);
  ck_assert_int_eq(ans, 0);
  ds.bits[3] = 0x001c0000;
  ck_assert(s21_from_decimal_to_int(ds, &ans) == 0);
  ck_assert_int_eq(ans, 0);
  ds.bits[2] = 0xffffffff;
  ds.bits[1] = 0xffffffff;
  ds.bits[0] = 0xffffffff;
  ds.bits[3] = 0x801c0000;

  ck_assert(s21_from_decimal_to_int(ds, &ans) == 0);
  ck_assert_int_eq(ans, -7);

  ds.bits[3] = 0x001c0000;
  ck_assert(s21_from_decimal_to_int(ds, &ans) == 0);
  ck_assert_int_eq(ans, 7);
}
END_TEST

START_TEST(CONV_TEST5) {
  s21_decimal ds = {0};
  ds.bits[0] = 0x7fffffff;
  int ans = 0;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, __INT32_MAX__);
  ds.bits[3] = 0x80000000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, -__INT32_MAX__);
  ds.bits[0] = 0x80000000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 1);
  ds.bits[3] = 0x00000000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 1);
  ds.bits[3] = 0x00010000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, __INT32_MAX__ / 10);
  ds.bits[3] = 0x80010000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, -__INT32_MAX__ / 10);
}
END_TEST

START_TEST(CONV_TEST6) {
  s21_decimal ds = {0};
  ds.bits[0] = 0x7fffffff;
  int ans = 0;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, __INT32_MAX__);
  ds.bits[3] = 0x80000000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, -__INT32_MAX__);
  ds.bits[0] = 0x80000000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 1);
  ds.bits[3] = 0x00000000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 1);
  ds.bits[3] = 0x00010000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, __INT32_MAX__ / 10);
  ds.bits[3] = 0x80010000;
  ck_assert_int_eq(s21_from_decimal_to_int(ds, &ans), 0);
  ck_assert_int_eq(ans, -__INT32_MAX__ / 10);
}
END_TEST

START_TEST(CONV_TEST7) {
  s21_decimal ds = {0};
  ds.bits[3] = 0x001c0000;
  float ans = 0.0;
  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  ck_assert(ans == 0.0);

  ds.bits[3] = 0x801c0000;
  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  ck_assert(ans == 0.0);

  ds.bits[3] = 0x00010000;
  ds.bits[0] = 0x00000040;
  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  // printf("%.50f\n", ans);
  ck_assert(ans == 6.4f);
  ds.bits[3] = 0x80010000;

  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  // printf("%.50f\n", ans);
  ck_assert(ans == -6.4f);

  ds.bits[3] = 0x00000000;
  ds.bits[0] = 0x25399253;
  ds.bits[1] = 0xaaaaaaaa;
  ds.bits[2] = 0x03523523;

  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  // printf("%.50f\n", ans);
  ck_assert(ans == 1027837890029898849765134931.0f);
  ds.bits[3] = 0x80000000;
  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  // printf("%.50f\n", ans);
  ck_assert(ans == -1027837890029898849765134931.0f);

  ds.bits[3] = 0x001c0000;
  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  // printf("%.50f\n", ans);
  ck_assert(ans == 0.1027837890029898849765134931f);
  ds.bits[3] = 0x801c0000;
  ck_assert_int_eq(s21_from_decimal_to_float(ds, &ans), 0);
  // printf("%.50f\n", ans);
  ck_assert(ans == -0.1027837890029898849765134931f);
}
END_TEST

Suite *dec_suite7(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal convert test");
  tc_core = tcase_create("FROM NUM TO DEC");
  tcase_add_test(tc_core, CONV_TEST1);
  tcase_add_test(tc_core, CONV_TEST2);
  tcase_add_test(tc_core, CONV_TEST3);
  tcase_add_test(tc_core, CONV_TEST4);
  tcase_add_test(tc_core, CONV_TEST5);
  tcase_add_test(tc_core, CONV_TEST6);
  tcase_add_test(tc_core, CONV_TEST7);
  suite_add_tcase(s, tc_core);

  return s;
}

///
// ROUND TESTS
/// \return

START_TEST(ROUND_TEST1) {
  s21_decimal ds = {0};
  s21_decimal res;
  ck_assert_int_eq(s21_truncate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ck_assert_int_eq(s21_round(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ck_assert_int_eq(s21_negate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ck_assert_int_eq(s21_floor(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ds.bits[3] == 0x801c0000;

  ck_assert_int_eq(s21_truncate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ck_assert_int_eq(s21_round(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ck_assert_int_eq(s21_negate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);

  ck_assert_int_eq(s21_floor(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0);
  ck_assert_int_eq(res.bits[3], 0);
}
END_TEST

START_TEST(ROUND_TEST2) {
  s21_decimal ds = {0};
  ds.bits[0] = 0x00000041;
  ds.bits[3] = 0x00010000;
  s21_decimal res;
  ck_assert_int_eq(s21_truncate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 6);
  ck_assert_int_eq(res.bits[3], 0);

  ck_assert_int_eq(s21_round(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 7);
  ck_assert_int_eq(res.bits[3], 0x00000000);

  ck_assert_int_eq(s21_negate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0x00000041);
  ck_assert_int_eq(res.bits[3], 0x80010000);

  ck_assert_int_eq(s21_floor(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0x00000006);
  ck_assert_int_eq(res.bits[3], 0);

  ds.bits[3] = 0x80010000;

  ck_assert_int_eq(s21_truncate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 6);
  ck_assert_int_eq(res.bits[3], 0x80000000);

  ck_assert_int_eq(s21_round(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 7);
  ck_assert_int_eq(res.bits[3], 0x80000000);

  ck_assert_int_eq(s21_negate(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 0x00000041);
  ck_assert_int_eq(res.bits[3], 0x00010000);

  ck_assert_int_eq(s21_floor(ds, &res), 0);
  ck_assert_int_eq(res.bits[2], 0);
  ck_assert_int_eq(res.bits[1], 0);
  ck_assert_int_eq(res.bits[0], 7);
  ck_assert_int_eq(res.bits[3], 0x80000000);
}
END_TEST

Suite *dec_suite8(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("Decimal ROUND test");
  tc_core = tcase_create("ROUND");
  tcase_add_test(tc_core, ROUND_TEST1);
  tcase_add_test(tc_core, ROUND_TEST2);
  suite_add_tcase(s, tc_core);

  return s;
}

int main() {
  int number_failed = EXIT_SUCCESS;
  Suite *s;
  SRunner *sr;

  ///
  s = dec_suite1();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);
  ///
  ///
  s = dec_suite2();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);
  ///
  ///
  s = dec_suite3();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);

  ///
  ///

  s = dec_suite4();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);

  ///
  ///

  s = dec_suite5();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);

  ///
  ///

  s = dec_suite6();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);

  s = dec_suite7();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);

  s = dec_suite8();
  sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = (number_failed == EXIT_FAILURE ? EXIT_FAILURE
                                                 : srunner_ntests_failed(sr));
  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
