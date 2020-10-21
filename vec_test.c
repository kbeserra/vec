#include "vec.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
VEC_DEF(1)

#include <stdlib.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#include <check.h>
#pragma clang diagnostic pop

START_TEST(projection_inbounds)
{
	vec1 a = { 1 };
	ck_assert_float_eq( vec1_pi( 0, a ), 1 );
	vec2 b = { 1, 2 };
	ck_assert_float_eq( vec2_pi( 0, b ), 1 );
	ck_assert_float_eq( vec2_pi( 1, b ), 2 );
	vec3 c = { 1, 2, 3 };
	ck_assert_float_eq( vec3_pi( 0, c ), 1 );
	ck_assert_float_eq( vec3_pi( 1, c ), 2 );
	ck_assert_float_eq( vec3_pi( 2, c ), 3 );
}
END_TEST


START_TEST(projection_outOfBounds)
{
	vec1 a = { 1 };
	ck_assert_float_nan( vec1_pi( 1, a ) );
	ck_assert_float_nan( vec1_pi( -1, a ) );
	vec2 b = { 1, 2 };
	ck_assert_float_nan( vec2_pi( 2, b ) );
	vec3 c = { 1, 2, 3 };
	ck_assert_float_nan( vec3_pi( 3, c ) );
}
END_TEST


START_TEST(add)
{
	vec1 a = { 1 };
	vec1_add( a, a, a);
	ck_assert( vec1_pi( 0, a ) == 2 );
	vec3 b = { 1, 2, 3 };
	vec3_add( b, b, b);
	ck_assert_float_eq( vec3_pi( 0, b ), 2 );
	ck_assert_float_eq( vec3_pi( 1, b ), 4 );
	ck_assert_float_eq( vec3_pi( 2, b ), 6 );
	vec3 c = { 6, 4, 2 };
	vec3 d;
	vec3_add( d, b, c);
	ck_assert_float_eq( vec3_pi( 0, d ), 8 );
	ck_assert_float_eq( vec3_pi( 1, d ), 8 );
	ck_assert_float_eq( vec3_pi( 2, d ), 8 );
}
END_TEST


START_TEST(sub)
{
	vec1 a = { 1 };
	vec1_sub( a, a, a);
	ck_assert_float_eq( vec1_pi( 0, a ), 0 );
	vec3 b = { 1, 2, 3 };
	vec3_sub( b, b, b);
	ck_assert_float_eq( vec3_pi( 0, b ), 0 );
	ck_assert_float_eq( vec3_pi( 1, b ), 0 );
	ck_assert_float_eq( vec3_pi( 2, b ), 0 );
	vec3 c = { 6, 4, 2 };
	vec3 d;
	vec3_sub( d, b, c);
	ck_assert_float_eq( vec3_pi( 0, d ), -6 );
	ck_assert_float_eq( vec3_pi( 1, d ), -4 );
	ck_assert_float_eq( vec3_pi( 2, d ), -2 );
}
END_TEST

START_TEST(scale)
{
	vec3 a = { 1, 2, 3 };
	vec3_scale( a, 1, a);
	ck_assert_float_eq( vec3_pi( 0, a ), 1 );
	ck_assert_float_eq( vec3_pi( 1, a ), 2 );
	ck_assert_float_eq( vec3_pi( 2, a ), 3 );
	vec3 b;
	vec3_scale( b, -1, a);
	ck_assert_float_eq( vec3_pi( 0, b ), -1 );
	ck_assert_float_eq( vec3_pi( 1, b ), -2 );
	ck_assert_float_eq( vec3_pi( 2, b ), -3 );
	vec3_scale( b, 100, a);
	ck_assert_float_eq( vec3_pi( 0, b ), 100 );
	ck_assert_float_eq( vec3_pi( 1, b ), 200 );
	ck_assert_float_eq( vec3_pi( 2, b ), 300 );
	vec3_scale( b, 1.0/2, a);
	ck_assert_float_eq( vec3_pi( 0, b ), 1.0/2 );
	ck_assert_float_eq( vec3_pi( 1, b ), 2.0/2 );
	ck_assert_float_eq( vec3_pi( 2, b ), 3.0/2 );
	vec3_scale( b, 1.0/3, a);
	ck_assert_float_eq( vec3_pi( 0, b ), 1.0/3 );
	ck_assert_float_eq( vec3_pi( 1, b ), 2.0/3 );
	ck_assert_float_eq( vec3_pi( 2, b ), 3.0/3 );
}
END_TEST

START_TEST(innerProduct)
{
	vec3 a = { 1, 2, 3 };
	ck_assert_float_eq( vec3_innerProduct( a, a), 14 );
	vec3 b;
	vec3_scale( b, -1, a);
	ck_assert_float_eq( vec3_innerProduct( a, b), -14 );
	ck_assert_float_eq( vec3_innerProduct( b, a), -14 );
	a[0] = 1; a[1] = a[2] = 0;
	b[1] = 1; b[0] = b[2] = 0;
	ck_assert_float_eq( vec3_innerProduct( b, a), 0 );
	ck_assert_float_eq( vec3_innerProduct( a, b), 0 );
}
END_TEST

START_TEST(norm)
{
	vec2 a = { cos(M_PI/4), sin(M_PI/4) };
	ck_assert_float_eq_tol( vec2_norm( a ), 1, VEC_EPSILON );
	vec2 b;
	vec2_scale( b, -1, a);
	ck_assert_float_eq_tol( vec2_norm( b ), 1, VEC_EPSILON );
	vec2_scale( b, -8, a);
	ck_assert_float_eq_tol( vec2_norm( b ), 8.0, VEC_EPSILON );
	vec2_scale( b, -16, a);
	ck_assert_float_eq_tol( vec2_norm( b ), 16.0, VEC_EPSILON );
}
END_TEST

START_TEST(normalize)
{
	vec2 a = { 0, 0 };
	vec2_normalize( a, a );
	ck_assert_float_nan( vec2_pi( 0, a ) );
	ck_assert_float_nan( vec2_pi( 1, a ) );

	srand(0);
	for( int i = 0; i < 1024; i++ ){
		a[0] = rand();
		a[1] = rand();
		vec2_normalize( a, a );
		if( fabs( a[0] ) < VEC_EPSILON && fabs( a[1]) < VEC_EPSILON ){
			ck_assert_float_nan( vec2_pi( 0, a ) );
			ck_assert_float_nan( vec2_pi( 1, a ) );
		} else {
			ck_assert_float_eq_tol( vec2_norm( a ), 1, VEC_EPSILON );
		}
	}
}
END_TEST




Suite * vec_test_suite(void){
	Suite *s;
	TCase *tc_core;

	s = suite_create("vec");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, projection_inbounds);
	tcase_add_test(tc_core, projection_outOfBounds);
	tcase_add_test(tc_core, add);
	tcase_add_test(tc_core, sub);
	tcase_add_test(tc_core, scale);
	tcase_add_test(tc_core, innerProduct);
	tcase_add_test(tc_core, norm);
	tcase_add_test(tc_core, normalize);
	suite_add_tcase(s, tc_core);

	return s;

}



MAT_DEF(3)


int vec_run_tests( ) {
	int number_failed;
	Suite *s;
	SRunner *sr;
	s = vec_test_suite();
	sr = srunner_create(s);
	srunner_run_all(sr, CK_NORMAL);
	number_failed = srunner_ntests_failed(sr);
	srunner_free(sr);

	return number_failed == 0;

}

#ifdef VEC_TEST

	int main( int argn, char** argv  ){
		return vec_run_tests() ? EXIT_SUCCESS : EXIT_FAILURE;
	}

#endif

#pragma clang diagnostic pop
