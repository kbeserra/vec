#ifndef VEC_H
#define VEC_H

#include <math.h>

#ifndef VEC_TYPE
  #define VEC_TYPE float
#endif

#ifndef VEC_FUNC
  #define VEC_FUNC static inline
#endif

#ifndef VEC_EPSILON
  #define VEC_EPSILON 1e-6
#endif

#ifndef VEC_SQRT
  #define VEC_SQRT sqrtl
#endif

#ifndef VEC_NAN
  #define VEC_NAN NAN
#endif


#define VEC_DEF(N)                                                             \
/* Type: vecN                                                                */\
/* ------------------                                                        */\
/* Member of VEC_TYPE^N viewed as a vector space. While the following is not */\
/* necessary, vecN is thought of as a column vector.                         */\
typedef VEC_TYPE vec##N[N];                                                    \
/* Function: vecN_pi                                                         */\
/* ------------------                                                        */\
/* Returns the projection of v onto the ith basis vector. Typically denoted  */\
/* as pi_i( v ) or v_i.                                                      */\
/*                                                                           */\
/* i: int, index of basis vector to project onto.                            */\
/* v: vecN, vector to project.                                               */\
/*                                                                           */\
/* returns: projection of v onto the ith basis vector of 0 <= i and i less   */\
/*      than the dimension of n. Otherwise, returns NAN.                     */\
VEC_FUNC VEC_TYPE vec##N##_pi( const int i, const vec##N v ) {                 \
  if( i < 0 || N <=  i ) {                                                     \
    return VEC_NAN;                                                            \
  }                                                                            \
  return (v[i]);                                                               \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_constant                                                   */\
/* ------------------                                                        */\
/* Sets each component of r to s.                                            */\
/*                                                                           */\
/* r: vecN, resulting vector.                                                */\
/* s: VEC_TYPE, scalar to populate r with.                                   */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_constant( vec##N r, const VEC_TYPE s ) {                \
  int i;                                                                       \
  for ( i = 0; i < N; i++ ) {                                                  \
    r[i] = s;                                                                  \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_zero                                                       */\
/* ------------------                                                        */\
/* Sets each component of r to zero.                                         */\
/*                                                                           */\
/* r: vecN, resulting zero vector.                                           */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_zero( vec##N r ) {                                      \
  vec##N##_constant( r, 0 );                                                   \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_nan                                                        */\
/* ------------------                                                        */\
/* Sets each component of r to nan, not a number.                            */\
/*                                                                           */\
/* r: vecN, resulting nan vector.                                            */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_nan( vec##N r ) {                                       \
  vec##N##_constant( r, VEC_NAN );                                             \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_lin3                                                       */\
/* ------------------                                                        */\
/* Computes the linear combination                                           */\
/*   r = q*u + s*v + t*w.                                                    */\
/*                                                                           */\
/* r: vecN, where the result is placed.                                      */\
/* q: VEC_TYPE, first scalar.                                                */\
/* u: vecN, first vector.                                                    */\
/* s: VEC_TYPE, second scalar.                                               */\
/* v: vecN, second vector.                                                   */\
/* t: VEC_TYPE, third scalar.                                                */\
/* w: vecN, third vector.                                                    */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_lin3( vec##N r, const VEC_TYPE q, const vec##N u,       \
                                       const VEC_TYPE s, const vec##N v,       \
                     const VEC_TYPE t, const vec##N w ) {                      \
  int i;                                                                       \
  for( i = 0; i < N; i++ ) {                                                   \
    r[i] = q * u[i] + s * v[i] + t * w[i];                                     \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_lin2                                                       */\
/* ------------------                                                        */\
/* Computes the linear combination                                           */\
/*   r = s*u + t*v.                                                          */\
/*                                                                           */\
/* r: vecN, where the result is placed.                                      */\
/* s: VEC_TYPE, first scalar.                                                */\
/* u: vecN, first vector.                                                    */\
/* t: VEC_TYPE, second scalar.                                               */\
/* v: vecN, second vector.                                                   */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_lin2( vec##N r, const VEC_TYPE s, const vec##N u,       \
                                       const VEC_TYPE t, const vec##N v ) {    \
  vec##N##_lin3(r, s, u, t, v, 0, u );                                         \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_add                                                        */\
/* ------------------                                                        */\
/* Computes the sum                                                          */\
/*   r = u + v.                                                              */\
/*                                                                           */\
/* r: vecN, where the result is placed.                                      */\
/* u: vecN, first vector.                                                    */\
/* v: vecN, second vector.                                                   */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_add( vec##N r, const vec##N u, const vec##N v ) {       \
  vec##N##_lin2( r, 1, u, 1, v );                                              \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_sub                                                        */\
/* ------------------                                                        */\
/* Computes the sum                                                          */\
/*   r = u - v.                                                              */\
/*                                                                           */\
/* r: vecN, where the result is placed.                                      */\
/* u: vecN, first vector.                                                    */\
/* v: vecN, second vector.                                                   */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_sub( vec##N r, const vec##N u, const vec##N v ) {       \
  vec##N##_lin2(r, 1, u, - 1, v);                                              \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_scale                                                      */\
/* ------------------                                                        */\
/* Computes the scalar product                                               */\
/*   r = s*v.                                                                */\
/*                                                                           */\
/* r: vecN, where the result is placed.                                      */\
/* s: VEC_TYPE, scalar.                                                      */\
/* v: vecN, vector.                                                          */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_scale( vec##N r, const VEC_TYPE s, const vec##N v) {    \
  vec##N##_lin2(r, s, v, 0, v);                                                \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_innerProduct                                               */\
/* ------------------                                                        */\
/* Computes the inner product typically denoted                              */\
/*    u cdot v = sum_i pi_i(u) * pi_i(v)                                     */\
/*             = sum_i u_i * v_i                                             */\
/*              = [ u, v ].                                                  */\
/*                                                                           */\
/* u: vecN, left hand vector in product.                                     */\
/* v: vecN, right hand vector in product.                                    */\
/*                                                                           */\
/* returns: VEC_TYPE, the result of the product.                             */\
VEC_FUNC VEC_TYPE vec##N##_innerProduct( const vec##N u, const vec##N v ) {    \
  VEC_TYPE rtn = 0;                                                            \
  int i;                                                                       \
  for( i = 0; i < N; i++ ) {                                                   \
    rtn += u[i] * v[i];                                                        \
  }                                                                            \
  return rtn;                                                                  \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_norm                                                       */\
/* ------------------                                                        */\
/* Computes the standard Euclidean norm of the vector v.                     */\
/*    | v | = sqrt( v cdot v )                                               */\
/*       = sqrt( sum_i pi_i(v) * pi_i(v) )                                   */\
/*       = sqrt( sum_i u_i * v_i ).                                          */\
/*                                                                           */\
/* v: vecN, vector to compute norm of.                                       */\
/*                                                                           */\
/* returns: VEC_TYPE, the Euclidean norm of v.                               */\
VEC_FUNC VEC_TYPE vec##N##_norm( const vec##N v) {                             \
  return ( (VEC_TYPE)VEC_SQRT( vec##N##_innerProduct(v, v) ) );                \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_metric                                                     */\
/* ------------------                                                        */\
/* Computes the standard Euclidean metric between the vectors u and v.       */\
/*    d(u,v) = | u - v |                                                     */\
/*                                                                           */\
/* u: vecN,                             .                                    */\
/* v: vecN,                             .                                    */\
/*                                                                           */\
/* returns: VEC_TYPE, the Euclidean distance between u and v.                */\
VEC_FUNC VEC_TYPE vec##N##_metric( const vec##N u, const vec##N v) {           \
  vec##N temp;                                                                 \
  vec##N##_sub( temp, u, v );                                                  \
  return vec##N##_norm( temp );                                                \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_normalize                                                  */\
/* ------------------                                                        */\
/* Computes the vector r with standard Euclidean norm 1 in the direction of  */\
/* the vector v.                                                             */\
/*    | r | = v / | v |.                                                     */\
/*                                                                           */\
/* r: vecN, unit vector in the direction of v.                               */\
/* v: vecN, vector to be normalized.                                         */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_normalize( vec##N r, const vec##N v ) {                 \
  VEC_TYPE n = vec##N##_norm(v);                                               \
  if( n < VEC_EPSILON ) {                                                      \
    vec##N##_nan( r );                                                         \
    return;                                                                    \
  }                                                                            \
  vec##N##_scale(r, 1.0 / n, v);                                               \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_project                                                    */\
/* ------------------                                                        */\
/* Computes the projection of u onto v and places the value into r           */\
/*    | r | = ( u dot v )/ ( v dot v ) v.                                    */\
/*                                                                           */\
/* r: vecN, vector in the direction v with norm the scalar projection of u   */\
/* onto v.                                                                   */\
/* u: vecN, vector being projected.                                          */\
/* v:  vecN, vector being projected onto.                                    */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_project( vec##N r, const vec##N u, const vec##N v ) {   \
  VEC_TYPE d = vec##N##_innerProduct(v, v);                                    \
  if( d < VEC_EPSILON ) {                                                      \
    vec##N##_nan( r );                                                         \
    return;                                                                    \
  }                                                                            \
  vec##N##_scale( r, vec##N##_innerProduct(u, v) / d, v );                     \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_max                                                        */\
/* ------------------                                                        */\
/* Computes the vector whose components are the maximum of the components of */\
/* u and v                                                                   */\
/*    r_i = max( u_i, v_i ).                                                 */\
/*                                                                           */\
/* r: vecN, vector whose components are the maximum of the components of u   */\
/* and v.                                                                    */\
/* u: vecN, vector.                                                          */\
/* v:  vecN, vector.                                                         */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_max( vec##N r, const vec##N u, const vec##N v ) {       \
  int i;                                                                       \
  for( i = 0; i < N; i++ ) {                                                   \
    r[i] = u[i] < v[i] ? v[i] : u[i];                                          \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \
/* Function: vecN_min                                                        */\
/* ------------------                                                        */\
/* Computes the vector whose components are the minimum of the components of */\
/* u and v                                                                   */\
/*    r_i = min( u_i, v_i ).                                                 */\
/*                                                                           */\
/* r: vecN, vector whose components are the minimum of the components of u   */\
/* and  v.                                                                   */\
/* u: vecN, vector.                                                          */\
/* v:  vecN, vector.                                                         */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void vec##N##_min( vec##N r, const vec##N u, const vec##N v ) {       \
  int i;                                                                       \
  for( i = 0; i < N; i++ ) {                                                   \
    r[i] = u[i] < v[i] ? u[i] : v[i];                                          \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \

VEC_DEF(2)

VEC_DEF(3)
// The following is a bit of a one-off.
// Function vec3_cross
// ------------------
// Computes the cross product of u and v
//   r = u \times v.
// The resulting vector r will be orthogonal to both u and v and have length
// equal to volume of the parallelogram spanned by the vectors u and v.
//
// r: vec3, the cross product of u and v.
// u: vec3, the left hand vector of the cross product.
// v: vec3, the right hand vector of the cross product.
VEC_FUNC void vec3_crossProduct( vec3 r, vec3 u, vec3 v ) {
  VEC_TYPE i, j, k;
  i = u[1] * v[2] - u[2] * v[1];
  j = u[2] * v[0] - u[0] * v[2];
  k = u[0] * v[1] - u[1] * v[0];
  r[0] = i;
  r[1] = j;
  r[2] = k;
}

VEC_DEF(4)


#define MAT_DEF(N)                                                             \
/* Type: matNxN                                                              */\
/* ------------------                                                        */\
/* matN is a linear transformation from vecN to vecN.                        */\
typedef vec##N mat##N##x##N[N];                                                \
/* Function: matN_identity                                                   */\
/* ------------------                                                        */\
/* Sets R to the identity transformation.                                    */\
/*                                                                           */\
/* R: matNxN, resulting transformation.                                      */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_identity( mat##N##x##N R ) {                      \
   int i, j;                                                                   \
   for( i = 0; i < N; i++ ) {                                                  \
     for( j = 0; j < N; j++ ) {                                                \
       R[i][j] = i == j ? 1.0 : 0.0;                                           \
     }                                                                         \
   }                                                                           \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_zero                                                     */\
/* ------------------                                                        */\
/* Sets R to the zero transformation.                                        */\
/*                                                                           */\
/* R: matNxN, resulting transformation.                                      */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_zero( mat##N##x##N R ) {                          \
   int i, j;                                                                   \
   for( i = 0; i < N; i++ ) {                                                  \
     for( j = 0; j < N; j++ ) {                                                \
       R[i][j] = 0.0;                                                          \
     }                                                                         \
   }                                                                           \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_nan                                                      */\
/* ------------------                                                        */\
/* Sets R to the nan transformation.                                         */\
/*                                                                           */\
/* R: matNxN, resulting transformation.                                      */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_nan( mat##N##x##N R ) {                           \
   int i;                                                                   \
   for( i = 0; i < N; i++ ) {                                                  \
     vec##N##_nan( R[i] );                                                     \
   }                                                                           \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_copy                                                     */\
/* ------------------                                                        */\
/* Copy the contents of M into R.                                            */\
/*                                                                           */\
/* R: matNxN, resulting transformation.                                      */\
/* M: matNxN, transformation to be coppied.                                  */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_copy( mat##N##x##N R, const mat##N##x##N M ) {    \
   int i, j;                                                                   \
   for( i = 0; i < N; i++ ) {                                                  \
     for( j = 0; j < N; j++ ) {                                                \
       R[i][j] = M[i][j];                                                      \
     }                                                                         \
   }                                                                           \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_row                                                      */\
/* ------------------                                                        */\
/* Obtains the ith row of the matrix M.                                      */\
/*                                                                           */\
/* r: vecN, desired row if matrix, if the passed index i is in [0,N),        */\
/*    otherwise r is set to the nan vector.                                  */\
/* M: matNxN, transformation to be coppied.                                  */\
/* i: int, index of row to obtain.                                           */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_row( vec##N r,                                    \
                                  const mat##N##x##N M, const int i ) {        \
 int j;                                                                        \
 if( i < 0 || N <= i ){                                                        \
   vec##N##_nan( r );                                                          \
   return  ;                                                                   \
 }                                                                             \
 for( j = 0; j<N; j++ ){                                                       \
   r[i] = M[j][i];                                                             \
 }                                                                             \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_column                                                   */\
/* ------------------                                                        */\
/* Obtains the ith column of the matrix M.                                   */\
/*                                                                           */\
/* r: vecN, desired column if matrix, if the passed index i is in [0,N),     */\
/*    otherwise r is set to the nan vector.                                  */\
/* M: matNxN, transformation to be coppied.                                  */\
/* i: int, index of column to obtain.                                        */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_column( vec##N r,                                 \
                                     const mat##N##x##N M, const int i ) {     \
 if( i < 0 || N <= i ){                                                        \
   vec##N##_nan( r );                                                          \
   return ;                                                                    \
 }                                                                             \
  vec##N##_scale( r, 1, M[i] );                                                \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_transpose                                                */\
/* ------------------                                                        */\
/* Computes the transpose of the matrix M.                                   */\
/*                                                                           */\
/* R: matNxN, the transpose of the matrix M.                                 */\
/* M: matNxN, matrix to compute transpose of.                                */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_transpose( mat##N##x##N R,                        \
                                        const mat##N##x##N M) {                \
  int i, j;                                                                    \
  VEC_TYPE tmp;                                                                \
  for( i = 0; i < N; i++ ){                                                    \
    for( j = 0; j <= i; j++ ){                                                 \
      tmp = M[i][j];                                                           \
      R[i][j] = M[j][i];                                                       \
      R[j][i] = tmp;                                                           \
    }                                                                          \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_lin2                                                     */\
/* ------------------                                                        */\
/* Computes linear combination                                               */\
/*   R = s * M + t * O                                                       */\
/*                                                                           */\
/* R: matNxN, result of the sum.                                             */\
/* s: VEC_TYPE, first scalar in the combination.                             */\
/* M: matNxN, first matrix in the combination.                               */\
/* t: VEC_TYPE, first scalar in the combination.                             */\
/* N: matNxN, second matrix in the combination.                              */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_lin2( mat##N##x##N R,                             \
                                   const VEC_TYPE s,                           \
                                   const mat##N##x##N M,                       \
                                   const VEC_TYPE t,                           \
                                   const mat##N##x##N O ) {                    \
  int i;                                                                       \
  for( i = 0; i < N; i++ ){                                                    \
    vec##N##_lin2( R[i], s, M[i], t, O[i] );                                   \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_add                                                      */\
/* ------------------                                                        */\
/* Computes sum                                                              */\
/*   R = M + O.                                                              */\
/*                                                                           */\
/* R: matNxN, result of the sum.                                             */\
/* M: matNxN, left summand.                                                  */\
/* O: matNxN, right summand.                                                 */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_add( mat##N##x##N R,                              \
                                  const mat##N##x##N M,                        \
                                  const mat##N##x##N O) {                      \
  mat##N##x##N##_lin2( R, 1, M, 1, O );                                        \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_sub                                                      */\
/* ------------------                                                        */\
/* Computes difference                                                       */\
/*   R = M - O.                                                              */\
/*                                                                           */\
/* R: matNxN, result of the difference.                                      */\
/* M: matNxN, left part of the sifference.                                   */\
/* O: matNxN, right part of the sifference.                                  */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_sub( mat##N##x##N R,                              \
                                  const mat##N##x##N M,                        \
                                  const mat##N##x##N O) {                      \
  mat##N##x##N##_lin2( R, 1, M, -1, O );                                       \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_scale                                                    */\
/* ------------------                                                        */\
/* Computes product                                                          */\
/*   R = s * M .                                                             */\
/*                                                                           */\
/* R: matNxN, result of the scale.                                           */\
/* s: VEC_TYPE, scalar.                                                      */\
/* M: matNxN, matrix.                                                        */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_scale( mat##N##x##N R,                            \
                                    const VEC_TYPE s,                          \
                                    const mat##N##x##N M ) {                   \
  mat##N##x##N##_lin2( R, s, M, 0, M );                                        \
}                                                                              \
                                                                               \
                                                                               \
/* Function: matNxN_multVec                                                  */\
/* ------------------                                                        */\
/* Computes the transformation of v my M                                     */\
/*   r = M v .                                                               */\
/* This function is not safe, if r and v share the same memory then the      */\
/* result, r, will be the zero vector.                                       */\
/*                                                                           */\
/* r: vecN, result of the transformation.                                    */\
/* M: matNxN, matrix.                                                        */\
/* v: vecN, vector to be transformed.                                        */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_multVec( vec##N r,                                \
                                      const mat##N##x##N M,                    \
                                      const vec##N v ) {                       \
  int i;                                                                       \
  vec##N##_zero( r );                                                          \
  for( i = 0; i < N; i++) {                                                    \
    vec##N##_lin2( r, 1, r, v[i], M[i] );                                      \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \
                                                                               \
                                                                               \
/* Function: matNxN_mult                                                     */\
/* ------------------                                                        */\
/* Computes the multiplaction my M * O                                       */\
/*   R = MO .                                                                */\
/* This function is not safe, each of the three matrices must have distinct  */\
/* memory.                                                                   */\
/*                                                                           */\
/* r: vecN, result of the transformation.                                    */\
/* M: matNxN, right matrix of the product.                                   */\
/* O: matNxN, left matrix of the product.                                    */\
/*                                                                           */\
/* returns: void                                                             */\
VEC_FUNC void mat##N##x##N##_mult( mat##N##x##N R,                             \
                                   const mat##N##x##N M,                       \
                                   const mat##N##x##N O ) {                    \
  int i;                                                                       \
  for( i = 0; i < N; i++) {                                                    \
    mat##N##x##N##_multVec( R[i], M, O[i] );                                   \
  }                                                                            \
}                                                                              \
                                                                               \
                                                                               \

MAT_DEF(2) 
MAT_DEF(3) 
MAT_DEF(4) 




#endif
