#ifndef S21_DECIMAL_H
#define S21_DECIMAL_H

#define FLOATSTRLEN 100
#define MAXFLOAT 10000000
#define MAXINT 9999999
#define HALFINT 16
#define SIGNOFS 31
#define ENDMASK 0x000000ff
#define MAXSCALE 28
#define BITSCNT 96
#define LASTBIT 0x80000000
#define LASRTHRBIT 0xe0000000
#define EVENPBITS 0xcccccccc
#define ODDPBITS 0x33333333
#define FIRSTWOBITS 0x00000003
#define UINT_MAX (1llu + __UINT32_MAX__)

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

typedef struct {
  uint32_t bits[4];
} s21_decimal;

void s21_inc(s21_decimal* dc);

int s21_from_float_to_decimal(float src, s21_decimal* dst);

int s21_get_sign(const s21_decimal dc);

void s21_change_sign(s21_decimal* dc);

int s21_from_int_to_decimal(int src, s21_decimal* dst);

int s21_set_scale(int scale, s21_decimal* dst);

void s21_set_zeros(s21_decimal* dst);

void s21_swap(s21_decimal* lhs, s21_decimal* rhs);

void s21_swap_int(int* lhs, int* rhs);

int s21_get_scale(const s21_decimal val);

int s21_try_upscale(s21_decimal* dc);

void s21_normalize(s21_decimal* dc);

int s21_add(s21_decimal value_1, s21_decimal value_2, s21_decimal* result);

int s21_sum(s21_decimal value1, s21_decimal value2, s21_decimal* result);

int s21_sub(s21_decimal value_1, s21_decimal value_2, s21_decimal* result);

void s21_normalize(s21_decimal* dc);

void s21_division2(s21_decimal* dc);

void s21_division5(s21_decimal* dc);

void s21_bank_round(s21_decimal* dc, int div);

void s21_bank_round_long(uint32_t* dc);

int s21_div2_long(uint32_t* dc);

void s21_inc_long(uint32_t* dc);

int s21_div5_long(uint32_t* dc);

void s21_division5_long(uint32_t* dc);

void s21_division2_long(uint32_t* dc);

int s21_div2(const s21_decimal dc);

int s21_div5(const s21_decimal dc);

int s21_is_less_abs(s21_decimal lhs, s21_decimal rhs);

int s21_is_less(s21_decimal lhs, s21_decimal rhs);

int s21_is_equal(s21_decimal lhs, s21_decimal rhs);

int s21_is_greater(s21_decimal lhs, s21_decimal rhs);

int s21_is_less_or_equal(s21_decimal lhs, s21_decimal rhs);

int s21_is_greater_or_equal(s21_decimal lhs, s21_decimal rhs);

int s21_is_not_equal(s21_decimal lhs, s21_decimal rhs);

void s21_move_point(s21_decimal* value_1, s21_decimal* value_2);

void s21_sub_abs(s21_decimal value_1, s21_decimal value_2, s21_decimal* result);

int s21_mul(s21_decimal value_1, s21_decimal value_2, s21_decimal* result);

void s21_fill_long(s21_decimal value_1, s21_decimal value_2, uint32_t* value);

int s21_div(s21_decimal value_1, s21_decimal value_2, s21_decimal* result);

int s21_from_decimal_to_int(s21_decimal src, int* dst);

int s21_from_decimal_to_float(s21_decimal src, float* dst);

int s21_floor(s21_decimal value, s21_decimal* result);

int s21_negate(s21_decimal value, s21_decimal* result);

int s21_truncate(s21_decimal value, s21_decimal* result);

int s21_round(s21_decimal value, s21_decimal* result);

void s21_bank_round_longl(unsigned long long* a, unsigned long long* b,
                          unsigned long long* c);

void s21_inc_longl(unsigned long long* a, unsigned long long* b,
                   unsigned long long* c);

void s21_reverse(s21_decimal* dc);

#endif