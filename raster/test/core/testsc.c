#include "rt_api.h"
#include "sc_raster.h"
#include "CUnit/Basic.h"
#include "cu_tester.h"

/* a test harness "includes" implementation which just
 * returns what it is told to return and logs the point
 * given to it.
 */
struct testset_inc {
	int includes ;
	LWPOINT *point ;
};

int testset_includes_fn(INCLUDES *inc, LWPOINT *p)
{
	struct testset_inc *ret ;

	CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(p) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(inc->params) ;

	ret = (struct testset_inc *)(inc->params) ;
	if (ret->point != NULL) {
		lwpoint_free(ret->point) ;
	}
	ret->point = lwgeom_as_lwpoint(lwgeom_clone_deep(lwpoint_as_lwgeom(p))) ;
	return ret->includes ;
}

INCLUDES *create_testset_includes(int retval)
{
	struct testset_inc *params ;
	INCLUDES *inc ;

	params = (struct testset_inc *)(lwalloc(sizeof (struct testset_inc)));
	CU_ASSERT_PTR_NOT_NULL_FATAL(params) ;

	params->includes = retval ;
	params->point = NULL ;

	inc = inc_create(params, testset_includes_fn, NULL) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;
	return inc ;
}

void destroy_testset_includes(INCLUDES *dead)
{
	if (!dead) {
		if (!dead->params) {
			struct testset_inc *p ;

			p = (struct testset_inc *)(dead->params) ;
			if (!p->point) {
				lwpoint_free(p->point) ;
			}

			lwfree(dead->params) ;
		}
		lwfree(dead) ;
	}
}

/*
 * a test harness "spatial collection" implementation which
 * binds the test-harness-includes with the mask-evaluator. This
 * gives the user control over the evaluated output.
 */
SPATIAL_COLLECTION *create_testset_sc(double inside, double outside)
{
	INCLUDES *inc ;
	EVALUATOR *eval ;
	SPATIAL_COLLECTION *sc ;
	GBOX fake_extent ; /* not used except as placeholder */

	/* instantiate a "testset" includes implementation */
	inc = create_testset_includes(1) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;

    /* instantiate the "mask evaluator". */
    eval = sc_create_mask_evaluator(inside, outside, 0) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(eval) ;

    /* instantiate a spatial collection from these */
    sc = sc_create(SPATIAL_PLUS_VALUE,SRID_UNKNOWN,&fake_extent,
    		NULL, inc, eval) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(sc) ;

    return sc ;
}

void destroy_testset_sc(SPATIAL_COLLECTION *dead)
{
	if (dead != NULL) {
		if (dead->inclusion != NULL) {
			destroy_testset_includes(dead->inclusion) ;
		}
		if (dead->evaluator != NULL) {
			sc_destroy_mask_evaluator(dead->evaluator) ;
		}
		sc_destroy(dead) ;
	}
}

void testset_includes_setval(INCLUDES *inc, int includes)
{
	struct testset_inc *params ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(inc->params) ;
	params = (struct testset_inc *)(inc->params) ;

	params->includes = includes ;
}

LWPOINT *testset_includes_getpoint(INCLUDES *inc)
{
	struct testset_inc *params ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(inc->params) ;
	params = (struct testset_inc *)(inc->params) ;

	return params->point ;
}

/*
 * An "EVALUATOR" implementation which records the point which was requested.
 * Used for testing the projection wrapper.
 */
struct testset_eval {
	LWPOINT *point ;
};

/*
 * An evaluator function which always returns the "result" as it is.
 * Prior to returning the result, it saves the point which the user
 * was evaluating.
 */
VALUE *testset_eval_fn(EVALUATOR *eval, LWPOINT *p)
{
	struct testset_eval *ret ;

	CU_ASSERT_PTR_NOT_NULL_FATAL(eval) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(p) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(eval->params) ;

	ret = (struct testset_eval *)(eval->params) ;
	if (ret->point != NULL) {
		lwpoint_free(ret->point) ;
	}
	ret->point = lwgeom_as_lwpoint(lwgeom_clone_deep(lwpoint_as_lwgeom(p))) ;
	return eval->result ;
}

EVALUATOR *create_testset_eval(void)
{
	struct testset_eval *params ;
	EVALUATOR *eval ;

	params = (struct testset_eval *)(lwalloc(sizeof (struct testset_eval)));
	CU_ASSERT_PTR_NOT_NULL_FATAL(params) ;

	params->point = NULL ;

	eval = eval_create(params, testset_eval_fn, NULL, 1) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(eval) ;
	return eval ;
}

void destroy_testset_eval(EVALUATOR *dead)
{
	if (!dead) {
		if (!dead->params) {
			struct testset_eval *p ;

			p = (struct testset_eval *)(dead->params) ;
			if (!p->point) {
				lwpoint_free(p->point) ;
			}

			lwfree(dead->params) ;
		}
		lwfree(dead) ;
	}
}

LWPOINT *testset_eval_getpoint(EVALUATOR *eval)
{
	struct testset_eval *params ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(eval) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(eval->params) ;
	params = (struct testset_eval *)(eval->params) ;

	return params->point ;
}

void test_dummy(void)
{
	printf("nothing to see here...\n");
}

/*
** The suite initialization function.
** Create any re-used objects.
*/
static int init_raster_sc_suite(void)
{
	return 0;
}

/*
** The suite cleanup function.
** Frees any global objects.
*/
static int clean_raster_sc_suite(void)
{
	return 0;
}


/*
** Used by test harness to register the tests in this file.
*/
CU_TestInfo raster_sc_tests[] =
{
	PG_TEST(test_dummy),
	CU_TEST_INFO_NULL
};
CU_SuiteInfo raster_sc_suite = {"Spatial Collection Test Suite (raster providers)",  init_raster_sc_suite,  clean_raster_sc_suite, raster_sc_tests};
