#include "f2c.h"
union {
    struct {
	doublereal a[676]	/* was [26][26] */, b[26], yx[26];
    } _1;
    struct {
	doublereal a[676]	/* was [26][26] */, b[26], y[26];
    } _2;
} lncoef_;
