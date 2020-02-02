/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"

#include "gmath.h"


/*******************************************************************************
 * Fixed size, vector and matrix operations useful for operations on 2D and 3D
 * coordinates.
 ******************************************************************************/


/*******************************************************************************
 * 2D vector operations
 ******************************************************************************/
void vec2d_zero(double a[2]) {

     a[0] = 0.0;
     a[1] = 0.0;
}



void vec2d_copy(double a2[2], double a1[2]) {

     a2[0] = a1[0];
     a2[1] = a1[1];
}



void vec2d_print(double a[2]) {

     printf("%.8f %.8f\n", a[0], a[1]);
}



void vec2d_scale(double s, double a[2], double b[2]) {

     b[0] = a[0] * s;
     b[1] = a[1] * s;
}



void vec2d_add(double a[2], double b[2], double c[2]) {

     c[0] = a[0] + b[0];
     c[1] = a[1] + b[1];
}



void vec2d_sub(double a[2], double b[2], double c[2]) {

     c[0] = a[0] - b[0];
     c[1] = a[1] - b[1];
}



void vec2d_lincmb(double a1[2], double r1, double a2[2], double r2, double b[2]) {

     b[0] = r1 * a1[0] + r2 * a2[0];
     b[1] = r1 * a1[1] + r2 * a2[1];
}



double vec2d_dot(double a[2], double b[2]) {

     return a[0]*b[0] + a[1]*b[1];
}



double vec2d_mag(double a[2]) {

     return sqrt(vec2d_dot(a, a));
}



void vec2d_unit(double a[2], double u[2]) {

     double mag;

     mag = vec2d_mag(a);

     u[0] = a[0] / mag;
     u[1] = a[1] / mag;
}



void vec2d_cross(double a[2], double b[2], double c[2]) {

}



void vec2d_dyadic(double a[2], double b[2], double c[2][2]) {

     c[0][0] = a[0] * b[0];
     c[0][1] = a[0] * b[1];

     c[1][0] = a[1] * b[0];
     c[1][1] = a[1] * b[1];
}



/*******************************************************************************
 * 3D vector operations
 ******************************************************************************/
void vec3d_zero(double a[3]) {

     a[0] = 0.0;
     a[1] = 0.0;
     a[2] = 0.0;
}



void vec3d_copy(double a2[3], double a1[3]) {

     a2[0] = a1[0];
     a2[1] = a1[1];
     a2[2] = a1[2];
}



void vec3d_print(double a[3]) {

     printf("%.8f %.8f %.8f\n", a[0], a[1], a[2]);
}



void vec3d_scale(double s, double a[3], double b[3]) {

     b[0] = a[0] * s;
     b[1] = a[1] * s;
     b[2] = a[2] * s;
}



void vec3d_add(double a[3], double b[3], double c[3]) {

     c[0] = a[0] + b[0];
     c[1] = a[1] + b[1];
     c[2] = a[2] + b[2];
}



void vec3d_sub(double a[3], double b[3], double c[3]) {

     c[0] = a[0] - b[0];
     c[1] = a[1] - b[1];
     c[2] = a[2] - b[2];
}



void vec3d_lincmb(double a1[3], double r1, double a2[3], double r2, double b[3]) {

     b[0] = r1 * a1[0] + r2 * a2[0];
     b[1] = r1 * a1[1] + r2 * a2[1];
     b[2] = r1 * a1[2] + r2 * a2[2];
}



double vec3d_dot(double a[3], double b[3]) {

     return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}



double vec3d_mag(double a[3]) {

     return sqrt(vec3d_dot(a, a));
}



void vec3d_unit(double a[3], double u[3]) {

     double mag;

     mag = vec3d_mag(a);

     u[0] = a[0] / mag;
     u[1] = a[1] / mag;
     u[2] = a[2] / mag;
}



void vec3d_cross(double a[3], double b[3], double c[3]) {

     c[0] = a[1] * b[2] - a[2] * b[1];
     c[1] = a[2] * b[0] - a[0] * b[2];
     c[2] = a[0] * b[1] - a[1] * b[0];
}



void vec3d_dyadic(double a[3], double b[3], double c[3][3]) {

     c[0][0] = a[0] * b[0];
     c[0][1] = a[0] * b[1];
     c[0][2] = a[0] * b[2];

     c[1][0] = a[1] * b[0];
     c[1][1] = a[1] * b[1];
     c[1][2] = a[1] * b[2];

     c[2][0] = a[2] * b[0];
     c[2][1] = a[2] * b[1];
     c[2][2] = a[2] * b[2];
}



/*******************************************************************************
 * 2D matrix operations
 ******************************************************************************/
void mat2d_zero(double a[2][2]) {

     a[0][0] = 0.0;
     a[0][1] = 0.0;

     a[1][0] = 0.0;
     a[1][1] = 0.0;
}



void mat2d_init(double a[2][2]) {

     a[0][0] = 1.0;
     a[0][1] = 0.0;

     a[1][0] = 0.0;
     a[1][1] = 1.0;
}



void mat2d_copy(double a2[2][2], double a1[2][2]) {

     a2[0][0] = a1[0][0];
     a2[0][1] = a1[0][1];

     a2[1][0] = a1[1][0];
     a2[1][1] = a1[1][1];
}



void mat2d_print(double a[2][2]) {

     printf("%.8f %.8f\n", a[0][0], a[0][1]);
     printf("%.8f %.8f\n", a[1][0], a[1][1]);
}



void mat2d_trans(double a[2][2], double b[2][2]) {

     b[0][0] = a[0][0];
     b[0][1] = a[1][0];

     b[1][0] = a[0][1];
     b[1][1] = a[1][1];
}



void mat2d_add(double a[2][2], double b[2][2], double c[2][2]) {

     c[0][0] = a[0][0] + b[0][0];
     c[0][1] = a[0][1] + b[0][1];

     c[1][0] = a[1][0] + b[1][0];
     c[1][1] = a[1][1] + b[1][1];
}



void mat2d_sub(double a[2][2], double b[2][2], double c[2][2]) {

     c[0][0] = a[0][0] - b[0][0];
     c[0][1] = a[0][1] - b[0][1];

     c[1][0] = a[1][0] - b[1][0];
     c[1][1] = a[1][1] - b[1][1];
}



void mat2d_scale(double s, double b[2][2], double c[2][2]) {

     c[0][0] = b[0][0] *s;
     c[0][1] = b[0][1] *s;

     c[1][0] = b[1][0] *s;
     c[1][1] = b[1][1] *s;
}



void mat2d_matvec(double a[2][2], double b[2], double c[2]) {

     c[0] = a[0][0]*b[0] + a[0][1]*b[1];
     c[1] = a[1][0]*b[0] + a[1][1]*b[1];
}



void mat2d_mul(double a[2][2], double b[2][2], double c[2][2]) {

     c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0];
     c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1];

     c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0];
     c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1];
}



int mat2d_inv(double a[2][2], double b[2][2]) {

     double det;

     det = mat2d_det(a);

     if (det != 0.0) {
          b[0][0] =  a[1][1] / det;
          b[0][1] = -a[0][1] / det;
          b[1][0] = -a[1][0] / det;
          b[0][0] =  a[0][0] / det;

          return 1;
     }
     else {
          return 0;
     }
}



double mat2d_det(double a[2][2]) {

     return a[0][0]*a[1][1] - a[0][1]*a[1][0];
}



void mat2d_rot(double a[2][2], double A) {

     double cosa;
     double sina;

     double temp_x1;

     cosa = cos(A);
     sina = sin(A);

     temp_x1 = a[0][0]*cosa - a[0][1]*sina;
     a[0][1] = a[0][0]*sina + a[0][1]*cosa;
     a[0][0] = temp_x1;

     temp_x1 = a[1][0]*cosa - a[1][1]*sina;
     a[1][1] = a[1][0]*sina + a[1][1]*cosa;
     a[1][0] = temp_x1;
}



/*******************************************************************************
 * 3D matrix operations
 ******************************************************************************/
void mat3d_zero(double a[3][3]) {

     a[0][0] = 0.0;
     a[0][1] = 0.0;
     a[0][2] = 0.0;

     a[1][0] = 0.0;
     a[1][1] = 0.0;
     a[1][2] = 0.0;

     a[2][0] = 0.0;
     a[2][1] = 0.0;
     a[2][2] = 0.0;
}



void mat3d_init(double a[3][3]) {

     a[0][0] = 1.0;
     a[0][1] = 0.0;
     a[0][2] = 0.0;

     a[1][0] = 0.0;
     a[1][1] = 1.0;
     a[1][2] = 0.0;

     a[2][0] = 0.0;
     a[2][1] = 0.0;
     a[2][2] = 1.0;
}



void mat3d_copy(double a2[3][3], double a1[3][3]) {

     a2[0][0] = a1[0][0];
     a2[0][1] = a1[0][1];
     a2[0][2] = a1[0][2];

     a2[1][0] = a1[1][0];
     a2[1][1] = a1[1][1];
     a2[1][2] = a1[1][2];

     a2[2][0] = a1[2][0];
     a2[2][1] = a1[2][1];
     a2[2][2] = a1[2][2];
}



void mat3d_print(double a[3][3]) {

     printf("%.8f %.8f %.8f\n", a[0][0], a[0][1], a[0][2]);
     printf("%.8f %.8f %.8f\n", a[1][0], a[1][1], a[1][2]);
     printf("%.8f %.8f %.8f\n", a[2][0], a[2][1], a[2][2]);
}



void mat3d_trans(double a[3][3], double b[3][3]) {

     b[0][0] = a[0][0];
     b[0][1] = a[1][0];
     b[0][2] = a[2][0];

     b[1][0] = a[0][1];
     b[1][1] = a[1][1];
     b[1][2] = a[2][1];

     b[2][0] = a[0][2];
     b[2][1] = a[1][2];
     b[2][2] = a[2][2];
}



void mat3d_add(double a[3][3], double b[3][3], double c[3][3]) {

     c[0][0] = a[0][0] + b[0][0];
     c[0][1] = a[0][1] + b[0][1];
     c[0][2] = a[0][2] + b[0][2];

     c[1][0] = a[1][0] + b[1][0];
     c[1][1] = a[1][1] + b[1][1];
     c[1][2] = a[1][2] + b[1][2];

     c[2][0] = a[2][0] + b[2][0];
     c[2][1] = a[2][1] + b[2][1];
     c[2][2] = a[2][2] + b[2][2];
}



void mat3d_sub(double a[3][3], double b[3][3], double c[3][3]) {

     c[0][0] = a[0][0] - b[0][0];
     c[0][1] = a[0][1] - b[0][1];
     c[0][2] = a[0][2] - b[0][2];

     c[1][0] = a[1][0] - b[1][0];
     c[1][1] = a[1][1] - b[1][1];
     c[1][2] = a[1][2] - b[1][2];

     c[2][0] = a[2][0] - b[2][0];
     c[2][1] = a[2][1] - b[2][1];
     c[2][2] = a[2][2] - b[2][2];
}



void mat3d_scale(double s, double b[3][3], double c[3][3]) {

     c[0][0] = b[0][0] *s;
     c[0][1] = b[0][1] *s;
     c[0][2] = b[0][2] *s;

     c[1][0] = b[1][0] *s;
     c[1][1] = b[1][1] *s;
     c[1][2] = b[1][2] *s;

     c[2][0] = b[2][0] *s;
     c[2][1] = b[2][1] *s;
     c[2][2] = b[2][2] *s;
}



void mat3d_matvec(double a[3][3], double b[3], double c[3]) {

     c[0] = a[0][0]*b[0] + a[0][1]*b[1] + a[0][2]*b[2];
     c[1] = a[1][0]*b[0] + a[1][1]*b[1] + a[1][2]*b[2];
     c[2] = a[2][0]*b[0] + a[2][1]*b[1] + a[2][2]*b[2];
}



void mat3d_mul(double a[3][3], double b[3][3], double c[3][3]) {

     c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
     c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
     c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];

     c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
     c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
     c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];

     c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
     c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
     c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
}



int mat3d_inv(double a[3][3], double b[3][3]) {

     double det;

     det = mat3d_det(a);

     if (det != 0.0) {
          b[0][0] = (a[1][1]*a[2][2] - a[1][2]*a[2][1]) / det;
          b[0][1] = (a[0][2]*a[2][1] - a[0][1]*a[2][2]) / det;
          b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1]) / det;

          b[1][0] = (a[1][2]*a[2][0] - a[1][0]*a[2][2]) / det;
          b[1][1] = (a[0][0]*a[2][2] - a[0][2]*a[2][0]) / det;
          b[1][2] = (a[0][2]*a[1][0] - a[0][0]*a[1][2]) / det;

          b[2][0] = (a[1][0]*a[2][1] - a[1][1]*a[2][0]) / det;
          b[2][1] = (a[0][1]*a[2][0] - a[0][0]*a[2][1]) / det;
          b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0]) / det;

          return 1;
     }
     else {
          return 0;
     }
}



double mat3d_det(double a[3][3]) {

     return  a[0][0]*a[1][1]*a[2][2] +
             a[0][1]*a[1][2]*a[2][0] +
             a[0][2]*a[1][0]*a[2][1] -
             a[0][2]*a[1][1]*a[2][0] -
             a[0][0]*a[1][2]*a[2][1] -
             a[0][1]*a[1][0]*a[2][2];
}



void mat3d_rot_x(double a[3][3], double x) {

     double cosx;
     double sinx;

     double temp_x2;

     cosx = cos(x);
     sinx = sin(x);

     temp_x2 = a[0][1]*cosx - a[0][2]*sinx;
     a[0][2] = a[0][1]*sinx + a[0][2]*cosx;
     a[0][1] = temp_x2;

     temp_x2 = a[1][1]*cosx - a[1][2]*sinx;
     a[1][2] = a[1][1]*sinx + a[1][2]*cosx;
     a[1][1] = temp_x2;

     temp_x2 = a[2][1]*cosx - a[2][2]*sinx;
     a[2][2] = a[2][1]*sinx + a[2][2]*cosx;
     a[2][1] = temp_x2;
}



void mat3d_rot_y(double a[3][3], double y) {

     double cosy;
     double siny;

     double temp_x1;

     cosy = cos(y);
     siny = sin(y);

     temp_x1 =  a[0][0]*cosy + a[0][2]*siny;
     a[0][2] = -a[0][0]*siny + a[0][2]*cosy;
     a[0][0] =  temp_x1;

     temp_x1 =  a[1][0]*cosy + a[1][2]*siny;
     a[1][2] = -a[1][0]*siny + a[1][2]*cosy;
     a[1][0] =  temp_x1;

     temp_x1 =  a[2][0]*cosy + a[2][2]*siny;
     a[2][2] = -a[2][0]*siny + a[2][2]*cosy;
     a[2][0] =  temp_x1;
}



void mat3d_rot_z(double a[3][3], double z) {

     double cosz;
     double sinz;

     double temp_x1;

     cosz = cos(z);
     sinz = sin(z);

     temp_x1 = a[0][0]*cosz - a[0][1]*sinz;
     a[0][1] = a[0][0]*sinz + a[0][1]*cosz;
     a[0][0] = temp_x1;

     temp_x1 = a[1][0]*cosz - a[1][1]*sinz;
     a[1][1] = a[1][0]*sinz + a[1][1]*cosz;
     a[1][0] = temp_x1;

     temp_x1 = a[2][0]*cosz - a[2][1]*sinz;
     a[2][1] = a[2][0]*sinz + a[2][1]*cosz;
     a[2][0] = temp_x1;
}



void mat3d_rot(double a[3][3], double x, double y, double z) {

     double cosx;
     double sinx;
     double cosy;
     double siny;
     double cosz;
     double sinz;

     double m11;
     double m12;
     double m13;
     double m21;
     double m22;
     double m23;
     double m31;
     double m32;
     double m33;

     double temp_x1;
     double temp_x2;

     cosx = cos(x);
     sinx = sin(x);
     cosy = cos(y);
     siny = sin(y);
     cosz = cos(z);
     sinz = sin(z);

     m11 =  cosy*cosz;
     m12 =  cosy*sinz;
     m13 = -siny;
     m21 =  sinx*siny*cosz - cosx*sinz;
     m22 =  sinx*siny*sinz + cosx*cosz;
     m23 =  sinx*cosy;
     m31 =  cosx*siny*cosz + sinx*sinz;
     m32 =  cosx*siny*sinz - sinx*cosz;
     m33 =  cosx*cosy;

     temp_x1 = a[0][0]*m11 + a[0][1]*m21 + a[0][2]*m31;
     temp_x2 = a[0][0]*m12 + a[0][1]*m22 + a[0][2]*m32;
     a[0][2] = a[0][0]*m13 + a[0][1]*m23 + a[0][2]*m33;
     a[0][0] = temp_x1;
     a[0][1] = temp_x2;

     temp_x1 = a[1][0]*m11 + a[1][1]*m21 + a[1][2]*m31;
     temp_x2 = a[1][0]*m12 + a[1][1]*m22 + a[1][2]*m32;
     a[1][2] = a[1][0]*m13 + a[1][1]*m23 + a[1][2]*m33;
     a[1][0] = temp_x1;
     a[1][1] = temp_x2;

     temp_x1 = a[2][0]*m11 + a[2][1]*m21 + a[2][2]*m31;
     temp_x2 = a[2][0]*m12 + a[2][1]*m22 + a[2][2]*m32;
     a[2][2] = a[2][0]*m13 + a[2][1]*m23 + a[2][2]*m33;
     a[2][0] = temp_x1;
     a[2][1] = temp_x2;
}
