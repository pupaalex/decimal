#include "s21_decimal.h"

int s21_from_float_to_decimal(float src, s21_decimal* dst) {
  s21_set_zeros(dst);
  int ret_code = 0;
  int scale = 0;
  if (src >= pow(2, BITSCNT) || src <= -pow(2, BITSCNT)) {
    ret_code = 1;
  }
  if ((src >= MAXFLOAT || src <= -MAXFLOAT) && ret_code == 0) {
    ret_code = s21_from_int_to_decimal((src > 0 ? MAXINT : -MAXINT), dst);
    s21_reverse(dst);
  } else if (ret_code == 0) {
    char buf[FLOATSTRLEN] = "";
    ret_code = (sprintf(buf, "%.40f", src) ? 0 : 1);

    int i = ((buf[0] == '+' || buf[0] == '-') ? 1 : 0);
    int flag = 0;
    for (; i < FLOATSTRLEN && buf[i] == '0'; ++i) {
    }
    flag = (buf[i] == '.');
    if (flag == 1 && ret_code == 0) {
      int zeros = 0;
      ++i;
      for (; i < FLOATSTRLEN && zeros < MAXSCALE && buf[i] == '0'; ++i) {
        ++zeros;
      }
      if (zeros > MAXSCALE - 1) {
        ret_code = 1;
      }
      scale = zeros;
    }
    if (ret_code == 0) {
      int a = 0;
      int cnt = 0;
      while (i < FLOATSTRLEN && cnt < 7) {
        if (isdigit(buf[i])) {
          a *= 10;
          a += (buf[i] - '0');
          ++cnt;
          if (flag == 1) {
            ++scale;
          }
        } else if (buf[i] == '.') {
          flag = 1;
        }
        ++i;
      }
      ret_code = s21_from_int_to_decimal((buf[0] == '-' ? -a : a), dst);
      s21_reverse(dst);
      while (scale > 0 && s21_div5(*dst) == 0 && s21_div2(*dst) == 0) {
        s21_division5(dst);
        s21_division2(dst);
        --scale;
      }
      ret_code = ((ret_code == 1) ? 1 : s21_set_scale(scale, dst));
    }
  }
  s21_normalize(dst);
  s21_reverse(dst);
  return ret_code;
}

int s21_get_sign(const s21_decimal dc) { return (dc.bits[3] >> SIGNOFS); }

void s21_change_sign(s21_decimal* dc) { dc->bits[3] ^= (1u << SIGNOFS); }

int s21_from_int_to_decimal(int src, s21_decimal* dst) {
  s21_set_zeros(dst);
  if (src < 0) {
    src = -src;
    s21_change_sign(dst);
  }
  dst->bits[2] = src;
  s21_reverse(dst);
  return 0;
}

int s21_set_scale(int scale, s21_decimal* dst) {
  int ret_code = 0;
  if (scale > MAXSCALE || scale < 0) {
    ret_code = 1;
  } else {
    dst->bits[3] &= (1u << SIGNOFS);
    dst->bits[3] += ((uint32_t)(scale << HALFINT));
  }
  return ret_code;
}

void s21_set_zeros(s21_decimal* dst) { memset(dst, 0, sizeof(int) * 4); }

int s21_get_scale(const s21_decimal val) {
  return (val.bits[3] >> HALFINT) & ENDMASK;
}

void s21_swap(s21_decimal* lhs, s21_decimal* rhs) {
  s21_decimal tmp = *rhs;
  *rhs = *lhs;
  *lhs = tmp;
}
void s21_swap_int(int* lhs, int* rhs) {
  int tmp = *rhs;
  *rhs = *lhs;
  *lhs = tmp;
}

int s21_try_upscale(s21_decimal* dc) {
  int ret_code = 0;
  s21_decimal tmp;
  s21_set_zeros(&tmp);
  if (((dc->bits[0] & LASRTHRBIT) != 0) || s21_get_scale(*dc) == 28) {
    ret_code = 1;
  } else {
    for (int i = 0; i < 3; ++i) {
      tmp.bits[i] = (dc->bits[i] << 3);
    }
    tmp.bits[0] += (dc->bits[1] >> 29);
    tmp.bits[1] += (dc->bits[2] >> 29);
    unsigned long long a = 0;
    for (int i = 2; i >= 0; --i) {
      a += ((unsigned long long)tmp.bits[i]) + (dc->bits[i] << 1);
      tmp.bits[i] = a % UINT_MAX;
      a /= UINT_MAX;
    }
    if (a > 0) {
      ret_code = 1;
    } else {
      a = ((unsigned long long)tmp.bits[1]) + (dc->bits[2] >> 31);
      tmp.bits[1] = a % UINT_MAX;
      a /= UINT_MAX;
      a += ((unsigned long long)tmp.bits[0]) + (dc->bits[1] >> 31);
      tmp.bits[0] = a % UINT_MAX;
      a /= UINT_MAX;
      if (a > 0) {
        ret_code = 1;
      } else {
        memcpy(dc->bits, tmp.bits, 3 * sizeof(uint32_t));
        int sc = s21_get_scale(*dc);
        s21_set_scale(sc + 1, dc);
      }
    }
  }
  return ret_code;
}

void s21_division2(s21_decimal* dc) {
  dc->bits[2] >>= 1;
  dc->bits[2] |= ((dc->bits[1] & 1) << 31);
  dc->bits[1] >>= 1;
  dc->bits[1] |= ((dc->bits[0] & 1) << 31);
  dc->bits[0] >>= 1;
}

void s21_normalize(s21_decimal* dc) {
  while (s21_div5(*dc) == 0 && s21_div2(*dc) == 0 && (s21_get_scale(*dc) > 0)) {
    s21_division2(dc);
    s21_division5(dc);
    s21_set_scale(s21_get_scale(*dc) - 1, dc);
  }
}

int s21_div5(const s21_decimal dc) {
  int sum = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 8; ++j) {
      sum += ((dc.bits[i] >> (j * 4)) & FIRSTWOBITS);
      sum -= ((dc.bits[i] >> (j * 4 + 2)) & FIRSTWOBITS);
    }
  }
  return (((sum % 5) + 5) % 5);
}

void s21_division5(s21_decimal* dc) {
  unsigned long long tmp = dc->bits[0] % 5;
  dc->bits[0] /= 5;
  tmp *= UINT_MAX;
  tmp += dc->bits[1];
  unsigned long long b = tmp % 5;
  dc->bits[1] = tmp / 5;
  b *= UINT_MAX;
  b += dc->bits[2];
  dc->bits[2] = b / 5;
}

int s21_div2(const s21_decimal dc) { return (dc.bits[2] % 2); }

void s21_inc(s21_decimal* dc) {
  if (dc->bits[2] == __UINT32_MAX__) {
    if (dc->bits[1] == __UINT32_MAX__) {
      dc->bits[0]++;
      dc->bits[1] = 0;
      dc->bits[2] = 0;
    } else {
      dc->bits[1]++;
      dc->bits[2] = 0;
    }
  } else {
    dc->bits[2]++;
  }
}

void s21_bank_round(s21_decimal* dc, int inv) {
  int flag = 0;
  if (s21_div2(*dc) != 0 || s21_div5(*dc) != 0) {
    if (s21_div5(*dc) == 0) {
      flag = 2;
    } else if (s21_div5(*dc) % 2 != s21_div2(*dc) && s21_div5(*dc) != 0) {
      flag = 1;
    }
  }
  s21_division2(dc);
  s21_division5(dc);
  if (flag == 1) {
    s21_inc(dc);
  } else if (flag == 2) {
    if (s21_div2(*dc) != inv) {
      s21_inc(dc);
    }
  }
}

void s21_bank_round_long(uint32_t* dc) {
  int flag = 0;
  if (s21_div2_long(dc) != 0 || s21_div5_long(dc) != 0) {
    if (s21_div5_long(dc) == 0) {
      flag = 2;
    } else if (s21_div5_long(dc) % 2 != s21_div2_long(dc) &&
               s21_div5_long(dc) != 0) {
      flag = 1;
    }
  }
  s21_division2_long(dc);
  s21_division5_long(dc);
  if (flag == 1) {
    s21_inc_long(dc);
  } else if (flag == 2) {
    if (s21_div2_long(dc) != 0) {
      s21_inc_long(dc);
    }
  }
}

void s21_bank_round_longl(unsigned long long* a, unsigned long long* b,
                          unsigned long long* c) {
  int ost = (*a) % 10;
  (*a) /= 10;

  (*b) += ost * UINT_MAX;

  ost = (*b) % 10;

  (*b) /= 10;

  (*c) += ost * UINT_MAX;
  ost = (*c) % 10;
  (*c) /= 10;

  if (ost > 5) {
    s21_inc_longl(a, b, c);
  } else if (ost == 5 && ((*c) % 2 == 1)) {
    s21_inc_longl(a, b, c);
  }
}

void s21_inc_longl(unsigned long long* a, unsigned long long* b,
                   unsigned long long* c) {
  if (*c == __UINT32_MAX__) {
    if (*b == __UINT32_MAX__) {
      (*a)++;
      (*b) = (*c) = 0;
    } else {
      (*b)++;
      (*c) = 0;
    }
  } else {
    (*c)++;
  }
}

int s21_div2_long(uint32_t* dc) { return (dc[5] % 2); }
int s21_div5_long(uint32_t* dc) {
  int sum = 0;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 8; ++j) {
      sum += ((dc[i] >> (j * 4)) & FIRSTWOBITS);
      sum -= ((dc[i] >> (j * 4 + 2)) & FIRSTWOBITS);
    }
  }
  return (((sum % 5) + 5) % 5);
}

void s21_division5_long(uint32_t* dc) {
  unsigned long long tmp = dc[0] % 5;
  dc[0] /= 5;
  tmp *= UINT_MAX;
  tmp += dc[1];
  unsigned long long b = tmp % 5;
  dc[1] = tmp / 5;
  b *= UINT_MAX;
  b += dc[2];
  unsigned long long c = b % 5;
  dc[2] = b / 5;
  c *= UINT_MAX;
  c += dc[3];
  unsigned long long d = c % 5;
  dc[3] = c / 5;
  d *= UINT_MAX;
  d += dc[4];
  unsigned long long e = d % 5;
  dc[4] = d / 5;
  e *= UINT_MAX;
  e += dc[5];
  dc[5] = e / 5;
}

void s21_division2_long(uint32_t* dc) {
  dc[5] >>= 1;
  dc[5] |= ((dc[4] & 1) << 31);
  dc[4] >>= 1;
  dc[4] |= ((dc[3] & 1) << 31);
  dc[3] >>= 1;
  dc[3] |= ((dc[2] & 1) << 31);
  dc[2] >>= 1;
  dc[2] |= ((dc[1] & 1) << 31);
  dc[1] >>= 1;
  dc[1] |= ((dc[0] & 1) << 31);
  dc[0] >>= 1;
}

void s21_inc_long(uint32_t* dc) {
  if (dc[5] == __UINT32_MAX__) {
    if (dc[4] == __UINT32_MAX__) {
      if (dc[3] == __UINT32_MAX__) {
        if (dc[2] == __UINT32_MAX__) {
          if (dc[1] == __UINT32_MAX__) {
            dc[0]++;
            dc[1] = dc[2] = dc[3] = dc[4] = dc[5] = 0;
          } else {
            dc[1]++;
            dc[2] = dc[3] = dc[4] = dc[5] = 0;
          }
        } else {
          dc[2]++;
          dc[3] = dc[4] = dc[5] = 0;
        }
      } else {
        dc[3]++;
        dc[4] = dc[5] = 0;
      }
    } else {
      dc[4]++;
      dc[5] = 0;
    }
  } else {
    dc[5]++;
  }
}

void s21_move_point(s21_decimal* value_1, s21_decimal* value_2) {
  int sc1 = s21_get_scale(*value_1);
  int sc2 = s21_get_scale(*value_2);
  int flag = 0;
  if (sc1 < sc2) {
    s21_swap(value_1, value_2);
    s21_swap_int(&sc1, &sc2);
    flag = 1;
  }

  while (sc1 > sc2) {
    int ret = s21_try_upscale(value_2);
    if (ret == 0) {
      sc2++;
    } else {
      break;
    }
  }
  while (sc1 > sc2 + 1) {
    s21_bank_round(value_1, 0);
    --sc1;
  }
  if (sc1 > sc2) {
    s21_bank_round(value_1, s21_div2(*value_2));
    --sc1;
  }
  s21_set_scale(sc1, value_1);
  s21_set_scale(sc2, value_2);
  if (flag) {
    s21_swap(value_1, value_2);
  }
}

int s21_add(s21_decimal value_1, s21_decimal value_2, s21_decimal* result) {
  int ret_code = 0;
  s21_reverse(&value_1);
  s21_reverse(&value_2);

  s21_normalize(&value_1);
  s21_normalize(&value_2);

  if (s21_get_sign(value_2) != s21_get_sign(value_1)) {
    if (s21_get_sign(value_1) == 0) {
      s21_change_sign(&value_2);
      s21_reverse(&value_1);
      s21_reverse(&value_2);
      ret_code = s21_sub(value_1, value_2, result);

    } else {
      s21_change_sign(&value_1);
      s21_reverse(&value_1);
      s21_reverse(&value_2);
      ret_code = s21_sub(value_2, value_1, result);
    }
    s21_reverse(result);
  } else {
    s21_set_zeros(result);
    s21_move_point(&value_1, &value_2);

    ret_code = s21_sum(value_1, value_2, result);
  }
  s21_normalize(result);
  s21_reverse(result);
  return ret_code;
}

int s21_sum(s21_decimal value_1, s21_decimal value_2, s21_decimal* result) {
  int ret_code = 0;

  unsigned long long a =
      ((unsigned long long)value_1.bits[2]) + value_2.bits[2];
  result->bits[2] = a % UINT_MAX;
  unsigned long long b =
      ((unsigned long long)value_1.bits[1]) + value_2.bits[1] + a / UINT_MAX;
  result->bits[1] = b % UINT_MAX;
  unsigned long long c =
      ((unsigned long long)value_1.bits[0]) + value_2.bits[0] + b / UINT_MAX;
  a = result->bits[2];
  b = result->bits[1];
  if (c > __UINT32_MAX__ && s21_get_scale(value_1) == 0) {
    if (s21_get_sign(value_1) == 0) {
      ret_code = 1;
    } else {
      ret_code = 2;
    }
  } else if (c > __UINT32_MAX__) {
    s21_bank_round_longl(&c, &b, &a);

    result->bits[0] = c;
    result->bits[1] = b;
    result->bits[2] = a;
    result->bits[3] = value_2.bits[3];
    s21_set_scale(s21_get_scale(value_1) - 1, result);

  } else {
    result->bits[0] = c;
    result->bits[3] = value_2.bits[3];
  }

  return ret_code;
}

int s21_sub(s21_decimal value_1, s21_decimal value_2, s21_decimal* result) {
  int ret_code = 0;
  s21_reverse(&value_1);
  s21_reverse(&value_2);
  s21_normalize(&value_1);
  s21_normalize(&value_2);

  if (s21_get_sign(value_2) != s21_get_sign(value_1)) {
    if (s21_get_sign(value_1) == 0) {
      s21_change_sign(&value_2);
      s21_reverse(&value_1);
      s21_reverse(&value_2);
      ret_code = s21_add(value_1, value_2, result);
    } else {
      s21_change_sign(&value_2);
      s21_reverse(&value_1);
      s21_reverse(&value_2);
      ret_code = s21_add(value_2, value_1, result);
    }
    s21_reverse(result);
  } else {
    uint32_t res_sign = s21_get_sign(value_1);
    if ((value_1.bits[3] >> 31) != 0) value_1.bits[3] -= (1u << 31);
    if ((value_2.bits[3] >> 31) != 0) value_2.bits[3] -= (1u << 31);
    s21_reverse(&value_1);
    s21_reverse(&value_2);
    if (s21_is_less(value_1, value_2)) {
      s21_swap(&value_1, &value_2);
      res_sign ^= 1;
    }
    s21_reverse(&value_1);
    s21_reverse(&value_2);

    s21_set_zeros(result);
    s21_move_point(&value_1, &value_2);

    s21_sub_abs(value_1, value_2, result);

    result->bits[3] |= (res_sign << SIGNOFS);
  }
  s21_normalize(result);
  s21_reverse(result);
  return ret_code;
}

void s21_sub_abs(s21_decimal value_1, s21_decimal value_2,
                 s21_decimal* result) {
  long long a = (unsigned long long)value_1.bits[2] - value_2.bits[2];

  int fg1 = (a < 0);
  if (fg1) {
    a += UINT_MAX;
  }
  result->bits[2] = a;
  long long b = (unsigned long long)value_1.bits[1] - value_2.bits[1] - fg1;
  int fg2 = (b < 0);
  if (fg2) {
    b += UINT_MAX;
  }
  result->bits[1] = b;
  result->bits[0] = value_1.bits[0] - value_2.bits[0] - fg2;
  s21_set_scale(s21_get_scale(value_2), result);
}

int s21_is_less_abs(s21_decimal lhs, s21_decimal rhs) {
  int sc1 = s21_get_scale(lhs);
  int sc2 = s21_get_scale(rhs);

  while (s21_try_upscale(&lhs) == 0) {
    ++sc1;
  }
  while (s21_try_upscale(&rhs) == 0) {
    ++sc2;
  }

  if (sc1 < sc2) {
    return 0;
  } else if (sc1 > sc2) {
    return 1;
  } else {
    return (lhs.bits[0] < rhs.bits[0]
                ? 1
                : (lhs.bits[0] > rhs.bits[0]
                       ? 0
                       : (lhs.bits[1] < rhs.bits[1]
                              ? 1
                              : (lhs.bits[1] > rhs.bits[1]
                                     ? 0
                                     : (lhs.bits[2] < rhs.bits[2] ? 1 : 0)))));
  }
}

int s21_is_less(s21_decimal lhs, s21_decimal rhs) {
  s21_reverse(&lhs);
  s21_reverse(&rhs);
  if (lhs.bits[0] == 0 && lhs.bits[1] == 0 && lhs.bits[2] == 0 &&
      rhs.bits[0] == 0 && rhs.bits[1] == 0 && rhs.bits[2] == 0) {
    return 0;
  }

  if (s21_get_sign(lhs) != s21_get_sign(rhs)) {
    if (s21_get_sign(lhs) == 0) {
      return 0;
    } else {
      return 1;
    }
  } else if (s21_get_sign(lhs) == 0) {
    return s21_is_less_abs(lhs, rhs);
  } else {
    return s21_is_less_abs(rhs, lhs);
  }
}

int s21_is_equal(s21_decimal lhs, s21_decimal rhs) {
  return ((!s21_is_less(lhs, rhs)) && (!s21_is_less(rhs, lhs)));
}

int s21_is_greater(s21_decimal lhs, s21_decimal rhs) {
  return ((!s21_is_less(lhs, rhs)) && (!s21_is_equal(rhs, lhs)));
}

int s21_is_less_or_equal(s21_decimal lhs, s21_decimal rhs) {
  return ((s21_is_less(lhs, rhs)) || (s21_is_equal(rhs, lhs)));
}

int s21_is_greater_or_equal(s21_decimal lhs, s21_decimal rhs) {
  return (!s21_is_less(lhs, rhs));
}

int s21_is_not_equal(s21_decimal lhs, s21_decimal rhs) {
  return s21_is_less(lhs, rhs) || s21_is_less(rhs, lhs);
}

int s21_mul(s21_decimal value_1, s21_decimal value_2, s21_decimal* result) {
  int ret_code = 0;
  s21_reverse(&value_1);
  s21_reverse(&value_2);
  s21_normalize(&value_1);
  s21_normalize(&value_2);

  uint32_t ans_sign = (s21_get_sign(value_1) ^ s21_get_sign(value_2));

  s21_set_zeros(result);
  int ans_sc = s21_get_scale(value_1) + s21_get_scale(value_2);
  uint32_t value[6];
  memset(value, 0, 6 * sizeof(uint32_t));

  s21_fill_long(value_1, value_2, value);

  while (value[0] != 0 || value[1] != 0 || value[2] != 0 || ans_sc > 28) {
    s21_bank_round_long(value);
    ans_sc--;
  }
  if (ans_sc < 0) {
    if (ans_sign == 0) {
      ret_code = 1;
    } else {
      ret_code = 2;
    }
  } else {
    memcpy(result->bits, value + 3, 3 * sizeof(uint32_t));
    s21_set_scale(ans_sc, result);
    result->bits[3] |= (ans_sign << SIGNOFS);
    s21_normalize(result);
  }
  s21_reverse(result);
  return ret_code;
}

void s21_fill_long(s21_decimal value_1, s21_decimal value_2, uint32_t* value) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      unsigned long long res =
          (unsigned long long)value_1.bits[i] * value_2.bits[j];
      for (int k = i + j + 1; k >= 0; --k) {
        unsigned long long ost = res / UINT_MAX;
        res %= UINT_MAX;
        res += value[k];
        value[k] = res % UINT_MAX;
        res /= UINT_MAX;
        res += ost;
      }
    }
  }
}

int s21_div(s21_decimal value_1, s21_decimal value_2, s21_decimal* result) {
  int ret_code = 0;
  s21_reverse(&value_1);
  s21_reverse(&value_2);
  s21_normalize(&value_1);
  s21_normalize(&value_2);
  s21_reverse(&value_1);
  s21_reverse(&value_2);
  s21_set_zeros(result);
  uint32_t res_sign = (s21_get_sign(value_1) ^ s21_get_sign(value_2));
  if (memcmp(value_2.bits, result->bits, 3 * sizeof(uint32_t)) == 0) {
    ret_code = 3;
  } else {
    if (s21_get_sign(value_1) == 1) value_1.bits[3] -= (1u << SIGNOFS);
    if (s21_get_sign(value_2) == 1) value_2.bits[3] -= (1u << SIGNOFS);
    s21_decimal tmp, ten, res, ans, tmp2, next, con, del;
    s21_from_int_to_decimal(0, &ans);
    s21_from_int_to_decimal(1, &tmp);
    s21_from_int_to_decimal(10, &ten);

    for (int i = 0; i < MAXSCALE; ++i) {
      s21_mul(tmp, ten, &res);
      memcpy(tmp.bits, res.bits, 4 * sizeof(uint32_t));
    }

    s21_from_float_to_decimal(0.1, &ten);
    int err = 0;
    do {
      if (s21_mul(tmp, value_2, &tmp2) == 0) {
        s21_set_zeros(&con);
        s21_set_zeros(&next);
        s21_set_zeros(&del);
        int r = s21_add(next, tmp2, &next);
        if (memcmp(tmp2.bits, result->bits, 3 * sizeof(uint32_t)) == 0) {
          r = 1;
          err = 1;
        }
        while ((r == 0) && s21_is_less_or_equal(next, value_1)) {
          memcpy(del.bits, next.bits, 4 * sizeof(uint32_t));
          r = s21_add(next, tmp2, &next);
          r = (s21_add(tmp, con, &con) ? 1 : r);
        }
        s21_sub(value_1, del, &res);
        memcpy(value_1.bits, res.bits, 4 * sizeof(uint32_t));
        ret_code = (s21_add(ans, con, &res) == 0 ? ret_code : 1);
        memcpy(ans.bits, res.bits, 4 * sizeof(uint32_t));
      }
      s21_mul(tmp, ten, &res);
      if (res.bits[0] == 0 && res.bits[1] == 0 && res.bits[2] == 0) {
        err = 1;
      } else {
        memcpy(tmp.bits, res.bits, 4 * sizeof(uint32_t));
      }
    } while (err == 0);

    ///
    s21_reverse(&tmp2);
    s21_try_upscale(&tmp2);
    s21_division2(&tmp2);
    s21_reverse(&tmp2);
    if (s21_is_less(tmp2, value_1)) {
      s21_add(ans, tmp, &ans);
    } else if (s21_is_equal(value_1, tmp2)) {
      if (ans.bits[0] % 2) {
        s21_add(ans, tmp, &ans);
      }
    }

    memcpy(result->bits, ans.bits, 4 * sizeof(uint32_t));
    result->bits[3] |= (res_sign << SIGNOFS);
    s21_reverse(result);
    s21_normalize(result);
  }
  s21_reverse(result);
  return ret_code;
}

int s21_from_decimal_to_int(s21_decimal src, int* dst) {
  int ret_code = 0;
  s21_reverse(&src);
  for (int i = 0; i < s21_get_scale(src); ++i) {
    s21_division2(&src);
    s21_division5(&src);
  }
  if (src.bits[0] != 0 || src.bits[1] != 0 || src.bits[2] > __INT32_MAX__) {
    ret_code = 1;
  } else {
    *dst = src.bits[2];
    if (s21_get_sign(src) == 1) {
      *dst = -(*dst);
    }
  }
  return ret_code;
}

int s21_from_decimal_to_float(s21_decimal src, float* dst) {
  s21_reverse(&src);
  char str[40];
  int ret_code = 0;
  memset(str, '0', 39);
  str[39] = '\0';
  int flag = 1;
  int sc = s21_get_scale(src);
  for (int i = 0; i < 29; ++i) {
    if (i == sc) {
      str[39 - i - flag] = '.';
      flag++;
    }
    int ost = 0;
    int f = s21_div5(src);
    int t = s21_div2(src);
    if (f % 2 != t % 2) ost += 5;
    ost += (f + 5) % 5;
    str[39 - i - flag] = (char)('0' + ost);
    s21_division2(&src);
    s21_division5(&src);
  }
  sscanf(str, "%f", dst);
  if (s21_get_sign(src) == 1) {
    (*dst) = -(*dst);
  }
  return ret_code;
}

int s21_floor(s21_decimal value, s21_decimal* result) {
  s21_reverse(&value);
  int ret_code = 0;
  if (value.bits[0] == 0 && value.bits[1] == 0 && value.bits[2] == 0) {
    memset(value.bits, 0, 4 * sizeof(uint32_t));
  }
  int sc = s21_get_scale(value);
  memcpy(result->bits, value.bits, 4 * sizeof(uint32_t));
  for (int i = 0; i < sc; ++i) {
    s21_division2(result);
    s21_division5(result);
  }
  s21_set_scale(0, result);
  s21_reverse(result);
  if (s21_get_sign(*result) == 1) {
    s21_decimal s;
    s21_from_int_to_decimal(1, &s);
    ret_code = s21_sub(*result, s, result);
  }
  return ret_code;
}

int s21_negate(s21_decimal value, s21_decimal* result) {
  int ret_code = 0;
  if (value.bits[0] == 0 && value.bits[1] == 0 && value.bits[2] == 0) {
    memset(value.bits, 0, 4 * sizeof(uint32_t));
  } else {
    memcpy(result->bits, value.bits, 4 * sizeof(uint32_t));
    s21_change_sign(result);
  }
  return ret_code;
}

int s21_truncate(s21_decimal value, s21_decimal* result) {
  int ret_code = 0;
  s21_reverse(&value);
  int sc = s21_get_scale(value);
  memcpy(result->bits, value.bits, 4 * sizeof(uint32_t));

  for (int i = 0; i < sc; ++i) {
    s21_division2(result);
    s21_division5(result);
  }
  s21_set_scale(0, result);
  s21_normalize(result);
  s21_reverse(result);
  return ret_code;
}

int s21_round(s21_decimal value, s21_decimal* result) {
  int ret_code = 0;
  s21_decimal tr, rn;
  s21_truncate(value, &tr);
  s21_sub(value, tr, &rn);
  s21_decimal half = {{0x00000005, 0x00000000, 0x00000000, 0x00010000}};
  s21_decimal neghalf = {{0x00000005, 0x00000000, 0x00000000, 0x80010000}};
  s21_decimal one = {{0x00000001, 0x00000000, 0x00000000, 0x00000000}};
  if (s21_is_greater_or_equal(rn, half)) {
    ret_code = s21_add(tr, one, result);
  } else if (s21_is_less_or_equal(rn, neghalf)) {
    ret_code = s21_sub(tr, one, result);
  } else {
    memcpy(result->bits, value.bits, 4 * sizeof(uint32_t));
  }
  return ret_code;
}

void s21_reverse(s21_decimal* dc) {
  uint32_t tmp = dc->bits[0];
  dc->bits[0] = dc->bits[2];
  dc->bits[2] = tmp;
}