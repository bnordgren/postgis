#include "rt_api.h"
#include "sc_raster.h"
#include "CUnit/Basic.h"
#include "cu_tester.h"
#include <math.h>
#include <stdio.h>

void
test_identity_gt(void)
{
	double o11, o12, o21, o22 ;
	double imag, jmag, theta_i, theta_ij ;

	/* first calculate the coefficients */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(1, 1, 0, M_PI_2, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, 1.0, 1e-6) ;

	/* then calculate the physically significant parameters. */
	imag = jmag = theta_i = theta_ij = -15.0 ;
	rt_raster_calc_phys_params(1, 0, 0, 1, &imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, M_PI_2, 1e-6) ;
}

void
test_bad_gt(void)
{
	double o11, o12, o21, o22 ;

	/* should fail on theta_ij == 180 degrees */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(1, 1, 0, M_PI, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, -15.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, -15.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, -15.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, -15.0, 1e-6) ;

	/* should fail on theta_ij == 0 degrees */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(1, 1, 0, 0, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, -15.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, -15.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, -15.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, -15.0, 1e-6) ;
}

void
test_flipped_j(void)
{
	double o11, o12, o21, o22 ;
	double imag, jmag, theta_i, theta_ij ;

	/* first calculate the coefficients */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(1, 1, 0, -M_PI_2, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, -1.0, 1e-6) ;

	/* then calculate the physically significant parameters. */
	imag = jmag = theta_i = theta_ij = -15.0 ;
	rt_raster_calc_phys_params(1, 0, 0, -1, &imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, -M_PI_2, 1e-6) ;
}

void
test_rotated(void)
{
	double o11, o12, o21, o22 ;
	double imag, jmag, theta_i, theta_ij ;

	/* first calculate the coefficients */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(1, 1, M_PI_4, M_PI_2, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, cos(M_PI_4), 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, sin(M_PI_4), 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, -sin(M_PI_4), 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, cos(M_PI_4), 1e-6) ;

	/* then calculate the physically significant parameters. */
	imag = jmag = theta_i = theta_ij = -15.0 ;
	rt_raster_calc_phys_params(cos(M_PI_4), sin(M_PI_4), -sin(M_PI_4), cos(M_PI_4),
			&imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 1.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, M_PI_4, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, M_PI_2, 1e-6) ;

}

void
test_scaled(void)
{
	double o11, o12, o21, o22 ;
	double imag, jmag, theta_i, theta_ij ;

	/* first calculate the coefficients */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(10, 20, 0, M_PI_2, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, 10.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, 20.0, 1e-6) ;

	/* then calculate the physically significant parameters. */
	imag = jmag = theta_i = theta_ij = -15.0 ;
	rt_raster_calc_phys_params(10, 0, 0, 20, &imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 10, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 20, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, M_PI_2, 1e-6) ;

}

void
test_scaled_flipped(void)
{
	double o11, o12, o21, o22 ;
	double imag, jmag, theta_i, theta_ij ;

	/* first calculate the coefficients */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(10, 20, 0, -M_PI_2, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, 10.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, -20.0, 1e-6) ;

	/* then calculate the physically significant parameters. */
	imag = jmag = theta_i = theta_ij = -15.0 ;
	rt_raster_calc_phys_params(10, 0, 0, -20, &imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 10, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 20, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, -M_PI_2, 1e-6) ;

}

void
test_skew(void)
{
	double o11, o12, o21, o22 ;
	double imag, jmag, theta_i, theta_ij ;
	double k_i, s_j ;

	k_i = tan(M_PI_4) ;
	s_j = 1./sqrt(k_i*k_i + 1) ;

	/* first calculate the coefficients */
	o11 = o12 = o21 = o22 = -15.0 ;
	rt_raster_calc_gt_coeff(1, 1, 0, M_PI_4, &o11, &o12, &o21, &o22) ;
	CU_ASSERT_DOUBLE_EQUAL(o11, 1, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o12, k_i*s_j, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o21, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(o22, s_j, 1e-6) ;

	/* then calculate the physically significant parameters. */
	imag = jmag = theta_i = theta_ij = -15.0 ;
	rt_raster_calc_phys_params(1, k_i*s_j, 0, s_j,
			&imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 1, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 1, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, M_PI_4, 1e-6) ;

}

void
test_raster_scale(void)
{
	rt_raster rast ;
	double imag, jmag, theta_i, theta_ij ;

	/* make the raster */
	rast = rt_raster_new(20,20) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(rast) ;

	/* first calculate the coefficients */
	rt_raster_set_phys_params(rast, 10, 20, 0, M_PI_2) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_x_scale(rast), 10.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_x_skew(rast), 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_y_skew(rast), 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_y_scale(rast), 20.0, 1e-6) ;

	/* then calculate the physically significant parameters. */
	rt_raster_get_phys_params(rast, &imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 10, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 20, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, M_PI_2, 1e-6) ;

	rt_raster_destroy(rast) ;
}

void
test_raster_skew(void)
{
	rt_raster rast ;
	double imag, jmag, theta_i, theta_ij ;
	double k_i, s_j ;

	/* precalculate terms */
	k_i = tan(M_PI_4) ;
	s_j = 1./sqrt(k_i*k_i + 1) ;

	/* make the raster */
	rast = rt_raster_new(10,10) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(rast) ;

	/* first calculate the coefficients */
	rt_raster_set_phys_params(rast, 1,1,0,M_PI_4) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_x_scale(rast), 1, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_x_skew(rast), k_i*s_j, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_y_skew(rast), 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(rt_raster_get_y_scale(rast), s_j, 1e-6) ;

	/* then calculate the physically significant parameters. */
	rt_raster_get_phys_params(rast, &imag, &jmag, &theta_i, &theta_ij) ;
	CU_ASSERT_DOUBLE_EQUAL(imag, 1, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(jmag, 1, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_i, 0.0, 1e-6) ;
	CU_ASSERT_DOUBLE_EQUAL(theta_ij, M_PI_4, 1e-6) ;

	rt_raster_destroy(rast);
}

/*
** The suite initialization function.
** Create any re-used objects.
*/
static int init_raster_gt_suite(void)
{
	return 0;
}

/*
** The suite cleanup function.
** Frees any global objects.
*/
static int clean_raster_gt_suite(void)
{
	return 0;
}


/*
** Used by test harness to register the tests in this file.
*/
CU_TestInfo raster_gt_tests[] =
{
	PG_TEST(test_identity_gt),
	PG_TEST(test_bad_gt),
	PG_TEST(test_flipped_j),
	PG_TEST(test_rotated),
	PG_TEST(test_scaled),
	PG_TEST(test_scaled_flipped),
	PG_TEST(test_skew),
	PG_TEST(test_raster_scale),
	PG_TEST(test_raster_skew),
	CU_TEST_INFO_NULL
};
CU_SuiteInfo raster_gt_suite = {"Raster Test Suite (geotransform calculations)",
		init_raster_gt_suite,  clean_raster_gt_suite, raster_gt_tests};
