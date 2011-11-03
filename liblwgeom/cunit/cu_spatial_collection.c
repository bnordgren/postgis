#include "CUnit/Basic.h"

#include "spatial_collection.h"
#include "liblwgeom_internal.h"
#include "cu_tester.h"
#include <math.h>



/*
 * tests for the correct behavior when converting a
 * "relationship string" to a code.
 */
static void test_get_relation_code(void)
{
	RELATION_TYPE reltype ;
	int success ;

	/* check bogus code */
	success = sc_get_relation_code("bogus", &reltype) ;
	CU_ASSERT_FALSE(success) ;

	/* check null */
	success = sc_get_relation_code("intersection", NULL) ;
	CU_ASSERT_FALSE(success) ;

	/* check intersection */
	success = sc_get_relation_code("intersection", &reltype) ;
	CU_ASSERT_TRUE(success) ;
	CU_ASSERT_EQUAL(reltype, INTERSECTION) ;

	/* check difference */
	success = sc_get_relation_code("difference", &reltype) ;
	CU_ASSERT_TRUE(success) ;
	CU_ASSERT_EQUAL(reltype, DIFFERENCE) ;

	/* check symmetric difference */
	success = sc_get_relation_code("symdifference", &reltype) ;
	CU_ASSERT_TRUE(success) ;
	CU_ASSERT_EQUAL(reltype, SYMDIFFERENCE) ;

	/* check union */
	success = sc_get_relation_code("union", &reltype) ;
	CU_ASSERT_TRUE(success) ;
	CU_ASSERT_EQUAL(reltype, UNION) ;

	/* uppercase fails */
	success = sc_get_relation_code("UNION", &reltype) ;
	CU_ASSERT_FALSE(success) ;
}


/*
 * direct references to the envelope code...
 */
extern GBOX *relation_env_intersection(LWPOLY *r1, LWPOLY *r2) ;
extern GBOX *relation_env_difference(LWPOLY *r1, LWPOLY *r2) ;
extern GBOX *relation_env_union(LWPOLY *r1, LWPOLY *r2) ;
extern GBOX *relation_env_symdifference(LWPOLY *r1, LWPOLY *r2) ;

/*
 * Tests that the correct envelope function is returned for
 * the given relationship code.
 */
void test_envelope_fn(void)
{
	ENVELOPE_PREP_OP env_fn ;

	env_fn = sc_get_envelope_fn(INTERSECTION) ;
	CU_ASSERT_EQUAL(env_fn, relation_env_intersection);

	env_fn = sc_get_envelope_fn(DIFFERENCE) ;
	CU_ASSERT_EQUAL(env_fn, relation_env_difference);

	env_fn = sc_get_envelope_fn(UNION) ;
	CU_ASSERT_EQUAL(env_fn, relation_env_union);

	env_fn = sc_get_envelope_fn(SYMDIFFERENCE) ;
	CU_ASSERT_EQUAL(env_fn, relation_env_symdifference);
}


/*
 * direct references to relationship function code
 */
extern int relation_intersection(int r1, int r2) ;
extern int relation_union(int r1, int r2) ;
extern int relation_difference(int r1, int r2) ;
extern int relation_symdifference(int r1, int r2) ;

/*
 * Tests that the correct "relation function" is returned for
 * the given relationship code.
 */
void test_relation_fn(void)
{
	RELATION_FN rel_fn ;

	rel_fn = sc_get_relation_fn(INTERSECTION) ;
	CU_ASSERT_EQUAL(rel_fn, relation_intersection) ;

	rel_fn = sc_get_relation_fn(UNION) ;
	CU_ASSERT_EQUAL(rel_fn, relation_union) ;

	rel_fn = sc_get_relation_fn(DIFFERENCE) ;
	CU_ASSERT_EQUAL(rel_fn, relation_difference) ;

	rel_fn = sc_get_relation_fn(SYMDIFFERENCE) ;
	CU_ASSERT_EQUAL(rel_fn, relation_symdifference) ;

}

/*
 * Sanity check on the behavior of the INTERSECTION relation
 * function
 */
void test_intersection_relation_fn(void)
{
	RELATION_FN rel_fn;

	rel_fn = sc_get_relation_fn(INTERSECTION) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(rel_fn) ;

	CU_ASSERT_FALSE(rel_fn(0,0)) ;
	CU_ASSERT_FALSE(rel_fn(0,1)) ;
	CU_ASSERT_FALSE(rel_fn(1,0)) ;
	CU_ASSERT_TRUE (rel_fn(1,1)) ;
}

/*
 * Sanity check on the behavior of the UNION relation
 * function
 */
void test_union_relation_fn(void)
{
	RELATION_FN rel_fn;

	rel_fn = sc_get_relation_fn(UNION) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(rel_fn) ;

	CU_ASSERT_FALSE(rel_fn(0,0)) ;
	CU_ASSERT_TRUE (rel_fn(0,1)) ;
	CU_ASSERT_TRUE (rel_fn(1,0)) ;
	CU_ASSERT_TRUE (rel_fn(1,1)) ;
}

/*
 * Sanity check on the behavior of the DIFFERENCE relation
 * function
 */
void test_difference_relation_fn(void)
{
	RELATION_FN rel_fn;

	rel_fn = sc_get_relation_fn(DIFFERENCE) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(rel_fn) ;

	CU_ASSERT_FALSE(rel_fn(0,0)) ;
	CU_ASSERT_FALSE(rel_fn(0,1)) ;
	CU_ASSERT_TRUE (rel_fn(1,0)) ;
	CU_ASSERT_FALSE(rel_fn(1,1)) ;
}

/*
 * Sanity check on the behavior of the SYMDIFFERENCE relation
 * function
 */
void test_symdifference_relation_fn(void)
{
	RELATION_FN rel_fn;

	rel_fn = sc_get_relation_fn(SYMDIFFERENCE) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(rel_fn) ;

	CU_ASSERT_FALSE(rel_fn(0,0)) ;
	CU_ASSERT_TRUE (rel_fn(0,1)) ;
	CU_ASSERT_TRUE (rel_fn(1,0)) ;
	CU_ASSERT_FALSE(rel_fn(1,1)) ;
}



/*
 * Support code to make a polygon in the shape of a box with
 * sides of the indicated length in the indicated SRID
 */
LWPOLY *
make_test_box(int srid, double length)
{
    POINTARRAY **rings = NULL;
    POINTARRAY *pts = NULL;
    LWPOLY* ret = NULL;
    POINT4D p4d;


    rings = (POINTARRAY **) lwalloc(sizeof (POINTARRAY*));
    CU_ASSERT_PTR_NOT_NULL_FATAL(rings) ;

    rings[0] = ptarray_construct(0, 0, 5);
    CU_ASSERT_PTR_NOT_NULL_FATAL(rings[0]) ;
    pts = rings[0];

    /* (0,0) first and last points */
    p4d.x = 0 ;
    p4d.y = 0 ;
    ptarray_set_point4d(pts, 0, &p4d);
    ptarray_set_point4d(pts, 4, &p4d); /* needed for closing it? */

    /* (0, length) */
    p4d.y = length ;
    ptarray_set_point4d(pts, 1, &p4d);

    /* (length, length) */
    p4d.x = length ;
    ptarray_set_point4d(pts, 2, &p4d);

    /* (length, 0) */
    p4d.y = 0 ;
    ptarray_set_point4d(pts, 3, &p4d);

    ret = lwpoly_construct(srid, 0, 1, rings);
    CU_ASSERT_PTR_NOT_NULL_FATAL(ret) ;

    return ret;
}


/*
 * Tests the includes implementation for a geometry wrapper
 * (e.g., inside the geometry is the "inside" value, outside the
 * geometry is the "outside" value)
 */
void test_geometry_wrapper_includes(void)
{
	LWGEOM *geometry ;
	INCLUDES *wrapper ;
	LWPOINT *test_pt ;
	int includes ;

	/* make a box and wrap it */
	geometry = lwpoly_as_lwgeom(make_test_box(SRID_UNKNOWN, 10.0)) ;
	wrapper = sc_create_geometry_includes(geometry) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(wrapper) ;

    /* index fn should be null */
    CU_ASSERT_PTR_NULL(wrapper->includesIndex) ;

    /* geopoint function should not be null */
    CU_ASSERT_PTR_NOT_NULL(wrapper->includes) ;

    /* test a point inside the box */
    test_pt = lwpoint_make2d(SRID_UNKNOWN, 5,5) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(test_pt) ;
    includes = wrapper->includes(wrapper, test_pt) ;
    CU_ASSERT_TRUE(includes) ;
    lwpoint_free(test_pt) ;

    /* test a point outside the box (x axis) */
    test_pt = lwpoint_make2d(SRID_UNKNOWN, 15,5) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(test_pt) ;
    includes = wrapper->includes(wrapper, test_pt) ;
    CU_ASSERT_FALSE(includes) ;
    lwpoint_free(test_pt) ;

    /* test a point outside the box (y axis) */
    test_pt = lwpoint_make2d(SRID_UNKNOWN, 5,15) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(test_pt) ;
    includes = wrapper->includes(wrapper, test_pt) ;
    CU_ASSERT_FALSE(includes) ;
    lwpoint_free(test_pt) ;

    /* test a point inside the box, except wrong srid */
    test_pt = lwpoint_make2d(24, 5,5) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(test_pt) ;
    includes = wrapper->includes(wrapper, test_pt) ;
    CU_ASSERT_FALSE(includes) ;
    lwpoint_free(test_pt) ;

	/* free memory */
    sc_destroy_geometry_includes(wrapper) ;
	lwgeom_free(geometry) ;
}

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
 * Tests the behavior of the "mask evaluator", which allows
 * the user to specify the value to return if the point is
 * inside (and a different value if the point is outside.)
 *
 * The mask evaluator is totally at the mercy of the "includes"
 * implementation of the spatial collection. If the includes
 * implementation claims that the point is included, the mask
 * evaluator should return the "inside" value. Otherwise it
 * should return the "outside" value. It does not attempt to
 * analyze the spatial relationship on its own.
 */
void test_mask_evaluator(void)
{
	SPATIAL_COLLECTION *sc ;
	LWPOINT *fake_point ;
	int includes ;
	VALUE *value ;
	double inside_val = 15 ;
	double outside_val = 20 ;

    /* instantiate a fake (ignored), non-null point */
    fake_point = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(fake_point) ;

    sc = create_testset_sc(inside_val, outside_val) ;

    /* test for "outside" value */
    testset_includes_setval(sc->inclusion, 0) ; /* set "includes" to false */
    includes = sc_includes(sc, fake_point) ;
    CU_ASSERT_FALSE(includes) ; /* verify not included */
    value = sc_evaluate(sc, fake_point) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(value) ;
    CU_ASSERT_EQUAL(value->length, 1) ;
    CU_ASSERT_EQUAL(value->data[0], outside_val) ;

    /* test for "inside" value */
    testset_includes_setval(sc->inclusion, 1) ; /* set "includes" to true */
    includes = sc_includes(sc, fake_point) ;
    CU_ASSERT_TRUE(includes) ; /* verify included */
    value = sc_evaluate(sc, fake_point) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(value) ;
    CU_ASSERT_EQUAL(value->length, 1) ;
    CU_ASSERT_EQUAL(value->data[0], inside_val) ;

    destroy_testset_sc(sc) ;
    lwpoint_free(fake_point) ;
}

/*
 * Tests the behavior of the "first value" evaluator. This evaluator
 * will return the value of the first spatial collection if it
 * exists, or the second one if not.
 *
 * The input values "exist" if the
 */
void test_first_value_evaluator(void)
{
	double first_inside = 1 ;
	double first_outside = 2 ;
	double second_inside = 3;
	double second_outside = 4 ;
	SPATIAL_COLLECTION *first ;
	SPATIAL_COLLECTION *second ;
	SPATIAL_COLLECTION *tester ;
	GBOX fake_extent ;
	EVALUATOR *eval ;
	INCLUDES *inc ;
	VALUE *value ;
	LWPOINT *fake_point ;

	/* instantiate a fake/ignored test point */
	fake_point = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(fake_point) ;

	/* create the two input collections */
	first = create_testset_sc(first_inside, first_outside) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(first) ;
	second = create_testset_sc(second_inside, second_outside);
    CU_ASSERT_PTR_NOT_NULL_FATAL(second) ;

    /* create the evaluator */
    eval = sc_create_first_value_evaluator(first, second) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(eval) ;

    /* create an includes instance */
    inc = create_testset_includes(1) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(inc) ;

    /* create the "testing" collection */
    tester = sc_twoinput_create(SPATIAL_PLUS_VALUE, NULL,
    		&fake_extent, inc, eval, first, second);
    CU_ASSERT_PTR_NOT_NULL_FATAL(tester) ;

    /* sanity checks */
    CU_ASSERT_TRUE(sc_hasTwoInputs(tester)) ;
    CU_ASSERT_PTR_NOT_NULL(tester->input1) ;
    CU_ASSERT_EQUAL(tester->input1, first) ;
    CU_ASSERT_PTR_NOT_NULL(tester->input2) ;
    CU_ASSERT_EQUAL(tester->input2, second) ;
    CU_ASSERT_PTR_NOT_NULL(tester->evaluator) ;
    CU_ASSERT_EQUAL(tester->evaluator, eval) ;
    CU_ASSERT_PTR_NOT_NULL(tester->inclusion) ;
    CU_ASSERT_EQUAL(tester->inclusion, inc) ;
    CU_ASSERT_PTR_NOT_NULL(eval->collection) ;
    CU_ASSERT_PTR_NOT_NULL(eval->result) ;
    CU_ASSERT_PTR_NOT_NULL(first->evaluator) ;
    CU_ASSERT_PTR_NOT_NULL(second->evaluator) ;
    CU_ASSERT_PTR_NOT_NULL(first->evaluator->result) ;
    CU_ASSERT_PTR_NOT_NULL(second->evaluator->result) ;
    CU_ASSERT_PTR_NOT_NULL(first->inclusion) ;
    CU_ASSERT_PTR_NOT_NULL(second->inclusion) ;
    CU_ASSERT_EQUAL(first->evaluator->result->length,
    		        second->evaluator->result->length) ;
    CU_ASSERT_EQUAL(eval->result->length,
    		        first->evaluator->result->length) ;



    /* neither input has a value */
    testset_includes_setval(first->inclusion, 0) ;
    testset_includes_setval(second->inclusion, 0) ;
    value = sc_evaluate(tester, fake_point) ;
    CU_ASSERT_PTR_NULL(value) ;

    /* first input has a value (only) */
    testset_includes_setval(first->inclusion, 1) ;
    CU_ASSERT_TRUE(sc_includes(first, fake_point)) ;
    CU_ASSERT_PTR_NOT_NULL(sc_evaluate(first, fake_point)) ;
    value = sc_evaluate(tester, fake_point) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(value) ;
    CU_ASSERT_EQUAL(value->length, 1) ;
    CU_ASSERT_EQUAL(value->data[0], first_inside) ;

    /* second input has a value (only) */
    testset_includes_setval(first->inclusion, 0) ;
    testset_includes_setval(second->inclusion, 1) ;
    value = sc_evaluate(tester, fake_point) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(value) ;
    CU_ASSERT_EQUAL(value->length, 1) ;
    CU_ASSERT_EQUAL(value->data[0], second_inside) ;

    /* both inputs have a value */
    testset_includes_setval(first->inclusion, 1) ;
    testset_includes_setval(second->inclusion, 1) ;
    value = sc_evaluate(tester, fake_point) ;
    CU_ASSERT_PTR_NOT_NULL_FATAL(value) ;
    CU_ASSERT_EQUAL(value->length, 1) ;
    CU_ASSERT_EQUAL(value->data[0], first_inside) ;

    /* clean up */
    sc_destroy_first_value_evaluator(eval) ;
    destroy_testset_sc(first) ;
    destroy_testset_sc(second) ;
    sc_twoinput_destroy(tester) ;
    destroy_testset_includes(inc) ;
    lwpoint_free(fake_point) ;
}


/*
 * some test harness infrastructure to ensure that the
 * "includes" implementation is retrieving the correct values...
 */

static int first_val = -1 ;
static int second_val = -1 ;

int testset_relation_op(int first, int second)
{
	first_val = first ;
	second_val = second ;
	return 0 ;  /* always return "not included"(false) */
}

/*
 * Tests the behavior of the "relation op" includes implementation.
 * This merely passes the "included" flag from the two inputs to
 * the "relation operator" The scope of this test is just to ensure
 * that these values are passed correctly.
 */
void test_relation_op_includes(void)
{
	INCLUDES *inc ;
	SPATIAL_COLLECTION *first ;
	SPATIAL_COLLECTION *second ;
	SPATIAL_COLLECTION *tester ;
	GBOX fake_extent ;
	LWPOINT *fake_point ;
	int included ;

	/* instantiate the fake test point */
	fake_point = lwpoint_make2d(SRID_UNKNOWN, 0, 0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(fake_point) ;

	/* create the two input collections. */
	first = create_testset_sc(1,2) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(first) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(first->inclusion) ;
	second = create_testset_sc(3,4) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(second) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(second->inclusion) ;

	/* create the relation op includes... */
	inc = sc_create_relation_includes(testset_relation_op) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(second) ;

	/* combine the inputs into a test set. */
	tester = sc_twoinput_create(SPATIAL_ONLY,
			NULL, &fake_extent, inc, NULL, first, second) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(tester) ;

	/* initialize the static testing variables */
	first_val = second_val = -1 ;

	/* neither input includes the test point */
	testset_includes_setval(first->inclusion, 0) ;
	testset_includes_setval(second->inclusion, 0) ;
	included = sc_includes(tester, fake_point) ;
	CU_ASSERT_EQUAL(included, 0) ; /* testset always returns false */
	CU_ASSERT_EQUAL(first_val, 0) ;
	CU_ASSERT_EQUAL(second_val, 0) ;

	/* first input (only) includes the test point */
	testset_includes_setval(first->inclusion, 1) ;
	testset_includes_setval(second->inclusion, 0) ;
	included = sc_includes(tester, fake_point) ;
	CU_ASSERT_EQUAL(first_val, 1) ;
	CU_ASSERT_EQUAL(second_val, 0) ;

	/* second input (only) includes the test point */
	testset_includes_setval(first->inclusion, 0) ;
	testset_includes_setval(second->inclusion, 1) ;
	included = sc_includes(tester, fake_point) ;
	CU_ASSERT_EQUAL(first_val, 0) ;
	CU_ASSERT_EQUAL(second_val, 1) ;

	/* both inputs include the test point */
	testset_includes_setval(first->inclusion, 1) ;
	testset_includes_setval(second->inclusion, 1) ;
	included = sc_includes(tester, fake_point) ;
	CU_ASSERT_EQUAL(first_val, 1) ;
	CU_ASSERT_EQUAL(second_val, 1) ;


	/* clean up */
	sc_twoinput_destroy(tester) ;
	destroy_testset_sc(first) ;
	destroy_testset_sc(second) ;
	lwpoint_free(fake_point) ;
}


#define SOURCE_P4TXT "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#define SOURCE_SRID 4326
#define DEST_P4TXT "+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#define DEST_SRID 32611

void test_includes_projection_wrapper(void)
{
	projPJ src ;
	projPJ dest ;
	LWPOINT *src_pt ;
	LWPOINT *dest_pt ;
	LWPOINT *test_pt ;
	POINT2D test2d, src2d ;
	INCLUDES *tester ;
	INCLUDES *wrapper ;
	int includes ;

	/* construct the projection objects... */
	src = lwproj_from_string(SOURCE_P4TXT) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(src) ;
	dest = lwproj_from_string(DEST_P4TXT) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(dest) ;

	/* Missoula, MT */
	src_pt = lwpoint_make2d(SOURCE_SRID,
			(-114.011593),
			46.862633);
	CU_ASSERT_PTR_NOT_NULL_FATAL(src_pt) ;

	dest_pt = lwgeom_as_lwpoint(lwgeom_clone_deep(lwpoint_as_lwgeom(src_pt))) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(dest_pt) ;
	lwgeom_transform(lwpoint_as_lwgeom(dest_pt), src, dest) ;
	dest_pt->srid = DEST_SRID ;

	/* make a test includes object to wrap */
	tester = create_testset_includes(0) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(tester) ;

	/* wrap it */
	wrapper = sc_create_projection_includes(tester, src, dest) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(wrapper) ;

	/* do we reproject? */
	/* since we wrapped a WGS84 collection and are presenting a
	 * utm zone 11 collection, we need to query in zone 11.
	 */
	testset_includes_setval(tester, 1) ;
	includes = wrapper->includes(wrapper, dest_pt) ;
	CU_ASSERT_TRUE(includes) ;

	/* however, the original collection has always been WGS84,
	 * so we expect the point presented to the original to have
	 * been projected.
	 */
	test_pt = testset_includes_getpoint(tester) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_pt) ;
	/* srid fields don't get changed when the points are projected,
	 * because Proj4 doesn't know anything about srids.
	 */
	CU_ASSERT_NOT_EQUAL(test_pt->srid, SOURCE_SRID) ;
	CU_ASSERT_EQUAL(test_pt->srid, DEST_SRID) ;

	/* so fix the srid before the compare... */
	test_pt->srid = SOURCE_SRID ;
	lwpoint_getPoint2d_p(test_pt, &test2d) ;
	lwpoint_getPoint2d_p(src_pt, &src2d) ;
	CU_ASSERT_DOUBLE_EQUAL(test2d.x, src2d.x, 0.0001) ;
	CU_ASSERT_DOUBLE_EQUAL(test2d.y, src2d.y, 0.0001) ;


	/* clean up */
	lwpoint_free(src_pt) ;
	lwpoint_free(dest_pt) ;
	pj_free(src) ;
	pj_free(dest) ;

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


void test_eval_projection_wrapper(void)
{
	projPJ src ;
	projPJ dest ;
	LWPOINT *src_pt ;
	LWPOINT *dest_pt ;
	LWPOINT *test_pt ;
	POINT2D test2d, src2d ;
	EVALUATOR *tester ;
	EVALUATOR *wrapper ;
	VALUE *val ;

	/* construct the projection objects... */
	src = lwproj_from_string(SOURCE_P4TXT) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(src) ;
	dest = lwproj_from_string(DEST_P4TXT) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(dest) ;

	/* Missoula, MT */
	src_pt = lwpoint_make2d(SOURCE_SRID,
			(-114.011593),
			46.862633);
	CU_ASSERT_PTR_NOT_NULL_FATAL(src_pt) ;

	dest_pt = lwgeom_as_lwpoint(lwgeom_clone_deep(lwpoint_as_lwgeom(src_pt))) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(dest_pt) ;
	lwgeom_transform(lwpoint_as_lwgeom(dest_pt), src, dest) ;
	dest_pt->srid = DEST_SRID ;

	/* make a test includes object to wrap */
	tester = create_testset_eval() ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(tester) ;

	/* wrap it */
	wrapper = sc_create_projection_eval(tester, src, dest) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(wrapper) ;

	/* do we reproject? */
	/* since we wrapped a WGS84 collection and are presenting a
	 * utm zone 11 collection, we need to query in zone 11.
	 */
	val = wrapper->evaluate(wrapper, dest_pt) ;
	CU_ASSERT_EQUAL(val, tester->result) ; /* return the wrapped result */

	/* however, the original collection has always been WGS84,
	 * so we expect the point presented to the original to have
	 * been projected.
	 */
	test_pt = testset_eval_getpoint(tester) ;
	CU_ASSERT_PTR_NOT_NULL_FATAL(test_pt) ;
	/* srid fields don't get changed when the points are projected,
	 * because Proj4 doesn't know anything about srids.
	 */
	CU_ASSERT_NOT_EQUAL(test_pt->srid, SOURCE_SRID) ;
	CU_ASSERT_EQUAL(test_pt->srid, DEST_SRID) ;

	/* so fix the srid before the compare... */
	test_pt->srid = SOURCE_SRID ;
	lwpoint_getPoint2d_p(test_pt, &test2d) ;
	lwpoint_getPoint2d_p(src_pt, &src2d) ;
	CU_ASSERT_DOUBLE_EQUAL(test2d.x, src2d.x, 0.0001) ;
	CU_ASSERT_DOUBLE_EQUAL(test2d.y, src2d.y, 0.0001) ;


	/* clean up */
	lwpoint_free(src_pt) ;
	lwpoint_free(dest_pt) ;
	pj_free(src) ;
	pj_free(dest) ;

}


/*
** The suite initialization function.
** Create any re-used objects.
*/
static int init_sc_suite(void)
{
	return 0;
}

/*
** The suite cleanup function.
** Frees any global objects.
*/
static int clean_sc_suite(void)
{
	return 0;
}

/*
** Used by test harness to register the tests in this file.
*/
CU_TestInfo sc_tests[] =
{
	PG_TEST(test_get_relation_code),
	PG_TEST(test_envelope_fn),
	PG_TEST(test_relation_fn),
	PG_TEST(test_geometry_wrapper_includes),
	PG_TEST(test_mask_evaluator),
	PG_TEST(test_first_value_evaluator),
	PG_TEST(test_relation_op_includes),
	PG_TEST(test_intersection_relation_fn),
	PG_TEST(test_union_relation_fn),
	PG_TEST(test_difference_relation_fn),
	PG_TEST(test_symdifference_relation_fn),
	PG_TEST(test_includes_projection_wrapper),
	PG_TEST(test_eval_projection_wrapper),
	CU_TEST_INFO_NULL
};
CU_SuiteInfo sc_suite = {"Spatial Collection Test Suite",  init_sc_suite,  clean_sc_suite, sc_tests};
