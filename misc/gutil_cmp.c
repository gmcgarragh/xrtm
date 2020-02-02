/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


#define cmp_frame(postfix_, type_)						\
int XCAT(cmp_, postfix_)(type_ a, type_ b) {					\
										\
     if (a < b)									\
          return -1;								\
     if (a > b)									\
          return 1;								\
										\
     return 0;									\
}


#define cmp_p_frame(postfix_, type_)						\
int XCAT(cmp_p_, postfix_)(type_ *a, type_ *b) {				\
										\
     if (*a < *b)								\
          return -1;								\
     if (*a > *b)								\
          return 1;								\
										\
     return 0;									\
}


#define cmp_n_frame(postfix_, type_)						\
int XCAT(cmp_n_, postfix_)(type_ *a, type_ *b, long n) {			\
										\
     long i;									\
										\
     for (i = 1; i < n; ++i) {							\
          if (a[i] < b[i])							\
               return -1;							\
          if (a[i] > b[i])							\
               return 1;							\
     }										\
										\
     return 0;									\
}


cmp_frame(uc, uchar)
cmp_n_frame(uc, uchar)
cmp_p_frame(uc, uchar)

cmp_frame(c, char)
cmp_n_frame(c, char)
cmp_p_frame(c, char)

cmp_frame(us, ushort)
cmp_n_frame(us, ushort)
cmp_p_frame(us, ushort)

cmp_frame(s, short)
cmp_n_frame(s, short)
cmp_p_frame(s, short)

cmp_frame(ui, uint)
cmp_n_frame(ui, uint)
cmp_p_frame(ui, uint)

cmp_frame(i, int)
cmp_n_frame(i, int)
cmp_p_frame(i, int)

cmp_frame(ul, ulong)
cmp_n_frame(ul, ulong)
cmp_p_frame(ul, ulong)

cmp_frame(l, long)
cmp_n_frame(l, long)
cmp_p_frame(l, long)

cmp_frame(ull, uint64_t)
cmp_n_frame(ull, uint64_t)
cmp_p_frame(ull, uint64_t)

cmp_frame(ll, int64_t)
cmp_n_frame(ll, int64_t)
cmp_p_frame(ll, int64_t)

cmp_frame(f, float)
cmp_n_frame(f, float)
cmp_p_frame(f, float)

cmp_frame(d, double)
cmp_n_frame(d, double)
cmp_p_frame(d, double)



#define cmp_dec_frame(postfix_, type_)						\
int XCAT(cmp_dec_, postfix_)(type_ a, type_ b) {					\
										\
     if (a > b)									\
          return -1;								\
     if (a < b)									\
          return 1;								\
										\
     return 0;									\
}


#define cmp_dec_p_frame(postfix_, type_)						\
int XCAT(cmp_dec_p_, postfix_)(type_ *a, type_ *b) {				\
										\
     if (*a > *b)								\
          return -1;								\
     if (*a < *b)								\
          return 1;								\
										\
     return 0;									\
}


#define cmp_dec_n_frame(postfix_, type_)						\
int XCAT(cmp_dec_n_, postfix_)(type_ *a, type_ *b, long n) {			\
										\
     long i;									\
										\
     for (i = 1; i < n; ++i) {							\
          if (a[i] > b[i])							\
               return -1;							\
          if (a[i] < b[i])							\
               return 1;							\
     }										\
										\
     return 0;									\
}


cmp_dec_frame(uc, uchar)
cmp_dec_n_frame(uc, uchar)
cmp_dec_p_frame(uc, uchar)

cmp_dec_frame(c, char)
cmp_dec_n_frame(c, char)
cmp_dec_p_frame(c, char)

cmp_dec_frame(us, ushort)
cmp_dec_n_frame(us, ushort)
cmp_dec_p_frame(us, ushort)

cmp_dec_frame(s, short)
cmp_dec_n_frame(s, short)
cmp_dec_p_frame(s, short)

cmp_dec_frame(ui, uint)
cmp_dec_n_frame(ui, uint)
cmp_dec_p_frame(ui, uint)

cmp_dec_frame(i, int)
cmp_dec_n_frame(i, int)
cmp_dec_p_frame(i, int)

cmp_dec_frame(ul, ulong)
cmp_dec_n_frame(ul, ulong)
cmp_dec_p_frame(ul, ulong)

cmp_dec_frame(l, long)
cmp_dec_n_frame(l, long)
cmp_dec_p_frame(l, long)

cmp_dec_frame(ull, uint64_t)
cmp_dec_n_frame(ull, uint64_t)
cmp_dec_p_frame(ull, uint64_t)

cmp_dec_frame(ll, int64_t)
cmp_dec_n_frame(ll, int64_t)
cmp_dec_p_frame(ll, int64_t)

cmp_dec_frame(f, float)
cmp_dec_n_frame(f, float)
cmp_dec_p_frame(f, float)

cmp_dec_frame(d, double)
cmp_dec_n_frame(d, double)
cmp_dec_p_frame(d, double)
