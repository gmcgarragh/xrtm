/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


void swap_data(void *a, void *b, void *buffer, size_t size) {

     memcpy(buffer, a, size);
     memcpy(a, b, size);
     memcpy(b, buffer, size);
}



void swap_index(int *a, int *b) {

    int temp;

    temp = *a;
    *a = *b;
    *b = temp;
}



void swap_pointer(uchar **a, uchar **b) {

    uchar *temp;

    temp = *a;
    *a = *b;
    *b = temp;
}



int insertion_sort_data(void *a, int n, size_t size,
                        int (*comp)(const void *, const void *)) {

     int j, p;

     void *temp;

     if (n > 1) {
          temp = (void *) malloc(size);
          if (temp == NULL) {
               fprintf(stderr, "ERROR: memory allocation failed\n");
               return -1;
          }

          for (p = 1; p < n; p++) {
               memcpy(temp, ((uchar *) a)+p*size, size);
               for (j = p; j > 0 && (*comp) (((uchar *) a)+(j*size-size), temp) > 0; j--) {
                    memcpy(((uchar *) a)+j*size, ((uchar *) a)+(j*size-size), size);
               }
               memcpy(((uchar *) a)+j*size, temp, size);
          }

          free(temp);
     }

     return 0;
}



void insertion_sort_index(void *a, int *index, int n, size_t size,
                          int (*comp)(const void *, const void *)) {

     int j, p, temp;

     for (p = 1; p < n; p++) {
          temp = index[p];
          for (j = p; j > 0 && (*comp) (((uchar *) a)+index[j-1]*size,
               ((uchar *) a)+temp*size) > 0; j--) {
               index[j] = index[j - 1];
          }
          index[j] = temp;
     }
}



void insertion_sort_pointer(void **a, int n, size_t size,
                        int (*comp)(const void *, const void *)) {

     int j, p;

     void *temp;

     for (p = 1; p < n; p++) {
          temp = *(((uchar **) a)+p);
          for (j = p; j > 0 && (*comp) (*(((uchar **) a)+(j-1)), temp) > 0; j--) {
               *(((uchar **) a)+j) = *(((uchar **) a)+(j-1));
          }
          *(((uchar **) a)+j) = (uchar *) temp;
     }
}



void *median_3_data(void *a, void *buffer, int lleft, int rright,
                    size_t size, int (*comp)(const void *, const void *)) {

     uchar *left, *right, *center;

     int ccenter = (lleft + rright) / 2;

     left = ((uchar *) a)+(size * lleft);
     right = ((uchar *) a)+(size * rright);
     center = ((uchar *) a)+(size * ccenter);

     if ((*comp) (left, center) > 0)
          swap_data(left, center, buffer, size);

     if ((*comp) (left, right) > 0)
          swap_data(left, right, buffer, size);

     if ((*comp) (center, right) > 0)
          swap_data(center, right, buffer, size);

     swap_data(center, right-size, buffer, size);

     return (void *) (right-size);
}



void q_sort_data(void *a, void *buffer, int left, int right,
                 size_t size, int (*comp)(const void *, const void *)) {

     int i, j;

     void *pivot;

     if (left + CUTOFF <= right) {
          pivot = median_3_data(a, buffer, left, right, size, comp);
          i = left;
          j = right - 1;

          for ( ; ; ) {
               while ((*comp) (((uchar *) a)+(++i*size), pivot) < 0) { }
               while ((*comp) (((uchar *) a)+(--j*size), pivot) > 0) { }
               if (i < j)
                   swap_data(((uchar *) a)+(i*size),
                             ((uchar *) a)+(j*size), buffer, size);
               else
                   break;
          }

          swap_data(((uchar *) a)+(i*size),
                   ((uchar *) a)+(right*size-size), buffer, size);

          q_sort_data(a, buffer, left, i - 1, size, comp);
          q_sort_data(a, buffer, i + 1, right, size, comp);
     }
     else
          insertion_sort_data(((uchar *) a)+(left*size),
                              right - left + 1, size, comp);
}



int quick_sort_data(void *a, int n, size_t size,
                    int (*comp)(const void *, const void *)) {

     void *buffer;

     if (n > CUTOFF) {
          buffer = (void *) malloc(size);
          if (buffer == NULL) {
               fprintf(stderr, "ERROR: memory allocation failed\n");
               return -1;
          }

          q_sort_data(a, buffer, 0, n - 1, size, comp);

          free(buffer);
     }
     else
          insertion_sort_data(a, n, size, comp);

     return 0;
}



void *median_3_index(void *a, int *index, int left, int right,
                     size_t size, int (*comp)(const void *, const void *)) {

     int center = (left + right) / 2;

     if ((*comp) (((uchar *) a)+(index[left]*size),
                  ((uchar *) a)+(index[center]*size)) > 0)
          swap_index(index + left, index + center);

     if ((*comp) (((uchar *) a)+(index[left]*size),
                  ((uchar *) a)+(index[right]*size)) > 0)
          swap_index(index + left, index + right);

     if ((*comp) (((uchar *) a)+(index[center]*size),
                  ((uchar *) a)+(index[right]*size)) > 0)
          swap_index(index + center, index + right);

     swap_index(index + center, index + right-1);

     return ((uchar *) a)+(index[right-1]*size);
}



void q_sort_index(void *a, int *index, int left, int right,
                  size_t size, int (*comp)(const void *, const void *)) {

     int i, j;

     void *pivot;

     if (left + CUTOFF <= right) {
          pivot = median_3_index(a, index, left, right, size, comp);
          i = left;
          j = right - 1;

          for ( ; ; ) {
               while ((*comp) (((uchar *) a)+(index[++i]*size), pivot) < 0) { }
               while ((*comp) (((uchar *) a)+(index[--j]*size), pivot) > 0) { }
               if (i < j)
                   swap_index(index + i, index + j);
               else
                   break;
          }

          swap_index(index + i, index + right-1);

          q_sort_index(a, index, left, i - 1, size, comp);
          q_sort_index(a, index, i + 1, right, size, comp);
     }
     else
          insertion_sort_index(a, index+left, right - left + 1, size, comp);
}



void quick_sort_index(void *a, int *index, int n, size_t size,
                      int (*comp)(const void *, const void *)) {

     if (n > CUTOFF)
          q_sort_index(a, index, 0, n - 1, size, comp);
     else
          insertion_sort_index(a, index, n, size, comp);
}



void **median_3_pointer(void **a, int lleft, int rright,
                    size_t size, int (*comp)(const void *, const void *)) {

     uchar **left, **right, **center;

     int ccenter = (lleft + rright) / 2;

     left = ((uchar **) a)+lleft;
     right = ((uchar **) a)+rright;
     center = ((uchar **) a)+ccenter;

     if ((*comp) (*left, *center) > 0)
          swap_pointer(left, center);

     if ((*comp) (*left, *right) > 0)
          swap_pointer(left, right);

     if ((*comp) (*center, *right) > 0)
          swap_pointer(center, right);

     swap_pointer(center, (right-1));

     return (void **) (right-1);
}



void q_sort_pointer(void **a, int left, int right,
                size_t size, int (*comp)(const void *, const void *)) {

     int i, j;

     void **pivot;

     if (left + CUTOFF <= right) {
          pivot = median_3_pointer(a, left, right, size, comp);
          i = left;
          j = right - 1;

          for ( ; ; ) {
               while ((*comp) (*(((uchar **) a)+(++i)), *pivot) < 0) { }
               while ((*comp) (*(((uchar **) a)+(--j)), *pivot) > 0) { }
               if (i < j)
                   swap_pointer(((uchar **) a)+i, ((uchar **) a)+j);
               else
                   break;
          }

          swap_pointer(((uchar **) a)+i, ((uchar **) a)+(right-1));

          q_sort_pointer(a, left, i - 1, size, comp);
          q_sort_pointer(a, i + 1, right, size, comp);
     }
     else
          insertion_sort_pointer((void **) (((uchar **) a)+left),
                                 right - left + 1, size, comp);
}



void quick_sort_pointer(void **a, int n, size_t size,
                    int (*comp)(const void *, const void *)) {

     if (n > CUTOFF)
          q_sort_pointer(a, 0, n - 1, size, comp);
     else
          insertion_sort_pointer(a, n, size, comp);
}
