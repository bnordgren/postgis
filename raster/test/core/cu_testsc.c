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

/*
 * a test harness "spatial collection" implementation which
 * binds the test-harness-includes with the test-harness-eval.
 */
SPATIAL_COLLECTION *create_testset_sc(void)
{
	INCLUDES *inc ;
	EVALUATOR *eval ;
	SPATIAL_COLLECTION *sc ;
	GBOX fake_extent ; /* not used except as placeholder */

	/* instantiate a "testset" includes implementation */
	inc = create_testset_includes(1) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;

    /* instantiate the "mask evaluator". */
    eval = create_testset_eval() ;
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
			destroy_testset_eval(dead->evaluator) ;
		}
		sc_destroy(dead) ;
	}
}

rt_raster make_test_raster(double width, double height,
		                   double xscale, double yscale,
		                   double xoff, double yoff)
{
	rt_raster grid_definition ;

	/* a 10x10 grid */
	grid_definition = rt_raster_new(width, height) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(grid_definition) ;

	/* scale 50x100 */
	rt_raster_set_scale(grid_definition, xscale, yscale) ;

	/* offset -10, +20 */
	rt_raster_set_offsets(grid_definition, xoff, yoff) ;

	return grid_definition ;
}

/*
 * Tests the ability to wrap a spatial collection such that it can be
 * accessed using the indices of a specified raster.
 */
void test_raster_aligned_collection(void)
{
	rt_raster grid_definition;
	SPATIAL_COLLECTION *tester ;
	SPATIAL_COLLECTION *aligned ;
	LWPOINT *test_point ;
	LWPOINT *probed_point ;

	/* 10x10 grid, 50x100 cell scale, -10,20 offset */
	grid_definition = make_test_raster(10,10, 50,100, -10,20) ;

	tester = create_testset_sc() ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(tester) ;
	aligned = sc_create_raster_aligned_collection(tester, grid_definition) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(aligned) ;

	/* test point should initially be null */
	test_point = testset_includes_getpoint(tester->inclusion) ;
	CU_ASSERT_PTR_NULL(test_point) ;
	test_point = testset_eval_getpoint(tester->evaluator) ;
	CU_ASSERT_PTR_NULL(test_point) ;

	/* ensure that we can access the collection normally (inclusion) */
	probed_point = lwpoint_make2d(SRID_UNKNOWN, 5,5) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(probed_point) ;
	sc_includes(aligned, probed_point) ;
	test_point = testset_includes_getpoint(tester->inclusion) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT(lwpoint_same(test_point, probed_point)) ; /* should be identical*/

	/* try accessing using indices defined by the grid (inclusion) */
	sc_includesIndex(aligned, probed_point) ;
	test_point = testset_includes_getpoint(tester->inclusion) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_DOUBLE_EQUAL(lwpoint_get_x(test_point),
			               (lwpoint_get_x(probed_point) * 50 - 10), 0.001) ;
	CU_ASSERT_DOUBLE_EQUAL(lwpoint_get_y(test_point),
			               (lwpoint_get_y(probed_point) * 100 +20), 0.001) ;

	/* ensure that we can access the collection normally (evaluator) */
	sc_evaluate(aligned, probed_point) ;
	test_point = testset_eval_getpoint(tester->evaluator) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT(lwpoint_same(test_point, probed_point)) ; /* should be identical*/

	/* try accessing using indices defined by the grid (evaluator) */
	sc_evaluateIndex(aligned, probed_point) ;
	test_point = testset_eval_getpoint(tester->evaluator) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_DOUBLE_EQUAL(lwpoint_get_x(test_point),
			               (lwpoint_get_x(probed_point) * 50 - 10), 0.001) ;
	CU_ASSERT_DOUBLE_EQUAL(lwpoint_get_y(test_point),
			               (lwpoint_get_y(probed_point) * 100 +20), 0.001) ;

	lwpoint_free(probed_point) ;
	sc_destroy_raster_aligned_collection(aligned) ;
	destroy_testset_sc(tester) ;
	rt_raster_destroy(grid_definition) ;
}

/*
 * Tests the "envelope" inclusion implementation. This one returns true
 * iff the point lies inside the grid defined by the raster. The comparison
 * is made in "index" space.
 */
void
test_envelope_inclusion(void)
{
	rt_raster grid ;
	INCLUDES *inc ;
	LWPOINT *test_point ;

	/* 50x50 grid; 20x40 cell scale; 60,-30 offset */
	grid = make_test_raster(50,50, 20,40, 60,-30) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(grid) ;

	/* make an "includes" which is based off of the envelope of the raster */
	inc = sc_create_raster_env_includes(grid) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;

	/* check that the "index" methods work appropriately...*/
	/* (0,0) is included */
	test_point = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_TRUE(inc->includesIndex(inc, test_point)) ;
	lwpoint_free(test_point) ;

	/* (-1,-1) is not included */
	test_point = lwpoint_make2d(SRID_UNKNOWN, -1,-1) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_FALSE(inc->includesIndex(inc, test_point)) ;
	lwpoint_free(test_point) ;

	/* (49,49) is included */
	test_point = lwpoint_make2d(SRID_UNKNOWN, 49,49) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_TRUE(inc->includesIndex(inc, test_point)) ;
	lwpoint_free(test_point) ;

	/* (50,50) is not included */
	test_point = lwpoint_make2d(SRID_UNKNOWN, 50,50) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_FALSE(inc->includesIndex(inc, test_point)) ;
	lwpoint_free(test_point) ;

	/* check that the "regular" methods (in the srid) work appropriately...*/
	/* (0,0) is included -> (60,-30) */
	test_point = lwpoint_make2d(SRID_UNKNOWN, 60,-30) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_TRUE(inc->includes(inc, test_point)) ;
	lwpoint_free(test_point) ;

	/* just short of (0,0) is not included */
	test_point = lwpoint_make2d(SRID_UNKNOWN, 59,-31) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_FALSE(inc->includes(inc, test_point)) ;
	lwpoint_free(test_point) ;

	/* Just short of (50,50) is included */
	test_point = lwpoint_make2d(SRID_UNKNOWN, 50*20+59,50*40-31) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_TRUE(inc->includes(inc, test_point)) ;
	lwpoint_free(test_point) ;

	/* (50,50) is not included */
	test_point = lwpoint_make2d(SRID_UNKNOWN, 50*20+60,50*40-30) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_point) ;
	CU_ASSERT_FALSE(inc->includes(inc, test_point)) ;
	lwpoint_free(test_point) ;

	sc_destroy_raster_env_includes(inc) ;
	rt_raster_destroy(grid) ;
}

/*
 * This adds the specified number of bands to the blank raster.
 * It initializes the band data with a very simple formula:
 *
 * (band, i, j) = band*width*height + j*width + i
 */
void
add_test_raster_data(rt_raster blank, int numbands)
{
	rt_band data ;
	int32_t *rawdata ;
	int width ;
	int height ;
	int i,j ;
	int band ;

	CU_ASSERT_PTR_NOT_NULL_FATAL(blank) ;
	width = rt_raster_get_width(blank) ;
	height = rt_raster_get_height(blank) ;

	for (band=0; band<numbands; band++){
		/* make a new data array */
		rawdata = rtalloc(sizeof(int32_t)*width*height) ;
		CU_ASSERT_PTR_NOT_NULL_FATAL(rawdata) ;

		/* initialize the array */
		for (j=0; j<height;j++) {
			for (i=0; i<width; i++) {
				rawdata[j*width+i] = band * width * height +
						             j * width  + i ;
			}
		}

		/* add the band data */
		data = rt_band_new_inline(width, height, PT_32BSI, 0,0, (uint8_t *)rawdata) ;
		CU_ASSERT_PTR_NOT_NULL_FATAL(rawdata) ;
		CU_ASSERT_NOT_EQUAL(rt_raster_add_band(blank,data,band), -1);
	}
}

/*
 * Test the ability to retrieve one or more band values per
 * location
 */
void
test_raster_bands_evaluator(void)
{
	rt_raster source ;
	EVALUATOR *eval ;
	VALUE *probe ;
	LWPOINT *pos ;
	int bands[] = {1,3} ;

	/* 5x5 grid, 1x1 cell scale, 0,0 offset */
	source = make_test_raster(5,5, 1,1, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(source) ;
	add_test_raster_data(source, 5) ; /* add 5 bands of data */

	/* check the "retrieve all bands" case */
	eval = sc_create_raster_bands_evaluator(source, NULL, 0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(eval) ;

	/* all five bands present? */
	CU_ASSERT_EQUAL(eval->result->length, 5) ;

	/* check to see that five bands are returned (index)*/
	pos = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(pos) ;
	probe = eval->evaluateIndex(eval, pos) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(probe) ;
	CU_ASSERT_EQUAL(probe->length, 5) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[0], 0, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[1], 25, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[2], 50, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[3], 75, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[4], 100, 0.01) ;
	lwpoint_free(pos) ;

	/* check to see that five bands are returned (geopoint)*/
	pos = lwpoint_make2d(SRID_UNKNOWN, 1,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(pos) ;
	probe = eval->evaluate(eval, pos) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(probe) ;
	CU_ASSERT_EQUAL(probe->length, 5) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[0], 1, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[1], 26, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[2], 51, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[3], 76, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[4], 101, 0.01) ;
	lwpoint_free(pos) ;

	sc_destroy_raster_bands_evaluator(eval) ;

	/* now select a subset of the available bands */
	eval = sc_create_raster_bands_evaluator(source, bands, 2) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(eval) ;

	/* should have 2 bands present */
	CU_ASSERT_EQUAL(eval->result->length, 2) ;

	/* check to see that two bands are returned (index) */
	pos = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(pos) ;
	probe = eval->evaluateIndex(eval, pos) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(probe) ;
	CU_ASSERT_EQUAL(probe->length, 2) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[0], 25, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[1], 75, 0.01) ;
	lwpoint_free(pos) ;

	/* check to see that two bands are returned (geopoint) */
	pos = lwpoint_make2d(SRID_UNKNOWN, 1,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(pos) ;
	probe = eval->evaluate(eval, pos) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(probe) ;
	CU_ASSERT_EQUAL(probe->length, 2) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[0], 26, 0.01) ;
	CU_ASSERT_DOUBLE_EQUAL(probe->data[1], 76, 0.01) ;
	lwpoint_free(pos) ;

	sc_destroy_raster_bands_evaluator(eval) ;
	rt_raster_destroy(source) ;
}


/*
 * checks to ensure that "includes" behaves correctly when the
 * raster contains a nodata value
 */
void
test_nodata_includes(void)
{
	rt_raster source ;
	rt_band   nodataband ;
	INCLUDES *nodata_inc ;
	LWPOINT *point ;

	/* 5x5 grid, 1x1 cell scale, 0,0 offset */
	source = make_test_raster(5,5, 1,1, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(source) ;
	add_test_raster_data(source, 5) ; /* add 5 bands of data */
	nodataband = rt_raster_get_band(source, 1) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(nodataband) ;
	rt_band_set_nodata(nodataband, 26) ; /* (1,0) is nodata */

	/* use band 1 as nodata reference */
	nodata_inc = sc_create_raster_nodata_includes(source, 1) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(nodata_inc) ;

	/* (0,0) is included (index) */
	point = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(point) ;
	CU_ASSERT_TRUE(nodata_inc->includesIndex(nodata_inc, point)) ;
	lwpoint_free(point) ;

	/* (0,0) is included (geopoint) */
	point = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(point) ;
	CU_ASSERT_TRUE(nodata_inc->includes(nodata_inc, point)) ;
	lwpoint_free(point) ;

	/* (1,0) is not included (index) */
	point = lwpoint_make2d(SRID_UNKNOWN, 1,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(point) ;
	CU_ASSERT_FALSE(nodata_inc->includesIndex(nodata_inc, point)) ;
	lwpoint_free(point) ;

	/* (1,0) is not included (geopoint) */
	point = lwpoint_make2d(SRID_UNKNOWN, 1,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(point) ;
	CU_ASSERT_FALSE(nodata_inc->includes(nodata_inc, point)) ;
	lwpoint_free(point) ;

	sc_destroy_raster_nodata_includes(nodata_inc) ;
	rt_raster_destroy(source) ;
}

/*
 * Tests to ensure that the sampling engine visits all the
 * pixels (where there is no nodata function), and produces the
 * correct result.
 */
void
test_sampling_engine_alldata(void)
{
	rt_raster source ;
	rt_raster result ;
	SPATIAL_COLLECTION *source_collection ;
	int i,j, band_k ; /* to loop over data */

	/* 5x5 grid, 1x1 cell scale, 0,0 offset */
	source = make_test_raster(5,5, 1,1, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(source) ;
	add_test_raster_data(source, 2) ; /* add 2 bands of data */
	source_collection = sc_create_raster_wrapper(source, 0, NULL, 0) ;

	/* make the result raster same dimensions and size as the above. */
	result = make_test_raster(5,5, 1,1, 0,0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(result) ;
	CU_ASSERT_NOT_EQUAL(rt_raster_get_num_bands(source),
			            rt_raster_get_num_bands(result)) ;

	/* run the sampling engine.
	 * this should copy all the bands and all the values to the result.
	 */
	sc_sampling_engine(source, result, NULL) ;

	CU_ASSERT_EQUAL(rt_raster_get_num_bands(source),
			        rt_raster_get_num_bands(result)) ;

	for (band_k=0; band_k<rt_raster_get_num_bands(source);band_k++) {
		rt_band source_band ;
		rt_band result_band ;

		/* get the band from the source */
		source_band = rt_raster_get_band(source, band_k) ;
		CU_ASSERT_PTR_NOT_NULL_FATAL(source_band) ;

		/* get the (hopefully copied) band from the result*/
		result_band = rt_raster_get_band(result, band_k) ;
		CU_ASSERT_PTR_NOT_NULL_FATAL(result_band) ;
		CU_ASSERT_NOT_EQUAL(source_band, result_band);


		for (j=0; j<rt_raster_get_height(source); j++) {
			for (i=0; i<rt_raster_get_width(source); i++) {
				double source_val ;
				double result_val ;

				rt_band_get_pixel(source_band,i,j,&source_val) ;
				rt_band_get_pixel(result_band,i,j,&result_val) ;
				CU_ASSERT_EQUAL(source_val, result_val) ;
			}
		}
	}

	rt_raster_destroy(source) ;
	rt_raster_destroy(result) ;
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
	PG_TEST(test_raster_aligned_collection),
	PG_TEST(test_envelope_inclusion),
	PG_TEST(test_nodata_includes),
	PG_TEST(test_raster_bands_evaluator),
	PG_TEST(test_sampling_engine_alldata),
	CU_TEST_INFO_NULL
};
CU_SuiteInfo raster_sc_suite = {"Spatial Collection Test Suite (raster providers)",  init_raster_sc_suite,  clean_raster_sc_suite, raster_sc_tests};
