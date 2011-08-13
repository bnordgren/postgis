#include "liblwgeom.h"
#include "spatial_collection.h"


/**
 * \defgroup includes_i Implementations of the "includes" interface.
 * @{
 */

/**
 * \defgroup sc_geom_include Inclusion in a single geometry
 * @{
 * Implementation of INCLUDES_FN which returns true if the evaluation
 * point is inside the collection's geometry.
 *
 */
int
geometry_includes(INCLUDES *inc, LWPOINT *point)
{
	LWGEOM *geom ; /* the geometry against which we are testing for inclusion*/

	if (inc == NULL || point == NULL) return 0 ;
	if (inc->params == NULL) return 0 ;

	geom = (LWGEOM *)(inc->params) ;

	return lwgeom_intersects(geom, lwpoint_as_lwgeom(point)) ;
}

/**
 * Constructor to create an instance of a "geometry_includes" implementation
 * of the Includes interface. This code keeps a reference to the "geom"
 * for evaluation, so do not change or destroy the "geom" until this
 * object is deallocated.
 * @param geom the geometry to evaluate.
 */
INCLUDES *
sc_create_geometry_includes(LWGEOM *geom)
{
	INCLUDES *gi ;

	gi = inc_create((PARAMETERS *)geom, geometry_includes, NULL) ;
	return gi ;
}


/**
 * Destroys an instance of the "geometry_includes" implementation
 * of the Includes interface. Only use this function to destroy objects
 * created with sc_create_geometry_includes.
 *
 */
void
sc_destroy_geometry_includes(INCLUDES *dead)
{
	inc_destroy(dead) ;
}

/** @} */ /* sc_geom_include documentation group */


/**
 * \defgroup sc_relation_include Inclusion for a combination of two collections.
 *
 * This group defines an implementation of an INCLUDES interface which is the
 * spatial combination of two input collections. The combination relationship
 * may be "union", "difference", "intersection", or "symdifference". This
 * INCLUDES interface operates on any combination of SPATIAL_COLLECTION types
 * (e.g., either input may be backed by a geometry, geometry-value, or a raster.)
 * @{
 */

/**
 * A RELATION_FN which evaluates the "intersection" of the two inputs.
 */
int relation_intersection(int r1, int r2)
{
	return r1 && r2 ;
}

/**
 * A RELATION_FN which evaluates the "union" of the two inputs.
 */
int relation_union(int r1, int r2)
{
	return r1 || r2 ;
}

/**
 * A RELATION_FN which evaluates the "difference" of the two inputs.
 */
int relation_difference(int r1, int r2)
{
	return r1 && !r2  ;
}

/**
 * A RELATION_FN which evaluates the symmetric difference of the
 * two inputs.
 */
int relation_symdifference(int r1, int r2)
{
	return (r1 && !r2) || (!r1 && r2) ;
}

/**
 * RELATION_FN is a function pointer typedef which describes a function
 * capable of evaluating a spatial relationship at a point. The arguments
 * are boolean values which reflect the value of the includes function
 * for the first and second inputs, respectively. The return value is
 * the result of the relationship evaluation.
 */
typedef int (RELATION_FN*)(int,int) ;


/**
 * Implementation of the INCLUDES_FN for use with a SPATIAL_COLLECTION
 * having two inputs. The INCLUDES_FN determines whether each of the
 * inputs includes the query point, then delegates the evaluation of
 * this collection's includes() value to the RELATION_FN specified
 * when the INCLUDES interface is created.
 */
int
relation_includes(INCLUDES *inc, LWPOINT *point)
{
	SPATIAL_COLLECTION *sc ;
	RELATION_FN relation ;
	int in1_included, in2_included ;

	if (inc == NULL || point == NULL || relation == NULL) return -1 ;
	if (inc->params == NULL) return -1 ;
	if (inc->collection == NULL) return -1 ;

	relation = (RELATION_FN)(inc->params) ;
	sc = inc->collection ;
	if (!sc_hasTwoInputs(sc)) return -1 ;

	in1_included = sc_includes(sc->input1, point) ;
	in2_included = sc_includes(sc->input2, point) ;

	return relation(in1_included, in2_included) ;
}

/**
 * Constructs an instance of INCLUDES which evaluates the combination
 * of two input collections.
 * @param relation the relation function used to evaluate the combination
 *        of the two input collections.
 */
INCLUDES *
sc_create_relation_includes(RELATION_FN relation)
{
	if (relation == NULL) return NULL ;
	return inc_create((PARAMETERS *)relation, relation_includes, NULL) ;
}

/**
 * Destroys an instance of the relation evaluator INCLUDES interface.
 * Only deallocate items created with sc_create_relation_includes.
 */
void
sc_destroy_relation_includes(INCLUDES *dead)
{
	inc_destroy(dead) ;
}

/** @}*/ /* end of sc_relation_include documentation group */

/** @}*/ /* end of includes_i documentation group */

/**
 * \defgroup evaluator_i Implementations of the "evaluator" interface.
 * @{
 */

/**
 * \defgroup mask_evaluator_i Evaluator interface returns a mask value
 *
 * This implementation examines the result of the collection's includes()
 * function. The user may set the numeric value associated with both
 * true and false values.
 * @{
 */

/**
 * Parametric data for the mask evaluator.
 */
struct mask_params_s {
	double false_val ;
	double true_val ;
};

/**
 * Workhorse function behind the mask evaluator. This function is
 * parameterized by an INCLUDE_FN in order to re-use the code for both
 * the evaluate() and evaluateIndex() functions.
 *
 * This function passes the query point to the specified includes()
 * function and returns a 1-band VALUE reflecting the result.
 */
VALUE *
mask_evaluator_util(EVALUATOR *eval, LWPOINT *point, INCLUDE_FN includes) {
	struct mask_params_s params ;
	VALUE *result ;
	double result_scalar ;

	/* check for naughty users */
	if (eval == NULL || point == NULL) return NULL  ;
	if (eval->collection == NULL) return NULL ;
	if (eval->result == NULL) return NULL ;
	if (eval->result->length != 1) return NULL ;

	params = (struct mask_params_s *)(eval->params) ;
	result = eval->result ;

	/* determine whether the collection includes the evaluation point or not */
	result_scalar = params->false_val ;
	if (includes(eval->collection->inclusion, point)) {
		result_scalar = params->true_val ;
	}
	result->data[0] = result_scalar ;

	return result ;
}

/**
 * Implementation of EVALUATOR_FN which leverages mask_evaluator_util
 * for the case of the evaluate() call. This requires that the utility
 * function be paired with the includes() method, as the point will be
 * in real-world coordinates and not raster indices.
 */
VALUE *
mask_evaluator(EVALUATOR *eval, LWPOINT *point)
{
	INCLUDE_FN includes ;

	/* check for naughty users */
	if (eval == NULL || point == NULL) return NULL  ;
	if (eval->collection == NULL) return NULL ;

	includes = eval->collection->inclusion->includes ;

	return mask_evaluator_util(eval, point, includes) ;

}

/**
 * Implementation of EVALUATOR_FN which leverages mask_evaluator_util
 * for the case of the evaluateIndex() call. This requires that the
 * utility function be paired with the includesIndex() method, as the
 * point will be expressed in raster indices instead of real-world coordinates.
 */
VALUE *
mask_index_evaluator(EVALUATOR *eval, LWPOINT *point)
{
	INCLUDE_FN includes ;

	/* check for naughty users */
	if (eval == NULL || point == NULL) return NULL  ;
	if (eval->collection == NULL) return NULL ;

	includes = eval->collection->inclusion->includesIndex ;
	if (includes == NULL ) return NULL ;

	return mask_evaluator_util(eval, point, includes) ;
}

/**
 * Constructor for the mask evaluation implementation of the Evaluator
 * interface.
 * @param true_val the value to return if this collection includes the
 *                 query point.
 * @param false_val the value to return if this collection does not
 *                  include the query point.
 * @param index flag controlling whether this instance will contain
 *              an evaluateIndex method.
 */
EVALUATOR *
sc_create_mask_evaluator(double true_val, double false_val, int index)
{
	struct mask_params_s *params ;
	EVALUATOR_FN idx_fn ;

	params = (struct mask_params_s *)lwalloc(sizeof(struct mask_params_s)) ;
	if (params == NULL) return NULL  ;

	params->true_val = true_val ;
	params->false_val = false_val ;

	/* decide whether or not to include the evaluateIndex() fn */
	idx_fn = NULL ;
	if (index) {
		idx_fn = mask_index_evaluator ;
	}

	return eval_create((PARAMETERS *)params, mask_evaluator, idx_fn) ;
}

/**
 * Destroys instances of mask evaluation objects. Only destroy objects
 * which were created with the sc_create_mask_evaluator function.
 */
void
sc_destroy_mask_evaluator(EVALUATOR *dead)
{
	if (dead == NULL) return ;
	if (dead->params != NULL) {
		lwfree(dead->params) ;
	}
	eval_destroy(dead) ;
}

/** @} */ /* end of mask_evaluator_i documentation group */


/**
 * \defgroup firstval_evaluator_i Evaluator to return a value from one of the collection's two inputs.
 *
 * This evaluator implements a simple default rule for determining the
 * value of a two-input collection. If the first input includes the query
 * point, then the value from the first input collection is returned.
 * Otherwise, the value from the second input collection is returned.
 *
 * This Evaluator requires that the two inputs have the same number of
 * bands. (e.g., the same "length" in their returned VALUEs).
 * @{
 */


/**
 * Utility function which
 */
VALUE *
first_value_util(EVALUATOR *eval,
		LWPOINT *point,
		INCLUDES_FN *inc1,
		INCLUDES_FN *inc2)
{
	VALUE *result, *result_in1, *result_in2 ;
	SPATIAL_COLLECTION *sc ;
	int have_value ;
	int i ;

	/* check for naughty users */
	if (eval == NULL || point == NULL || inc1 == NULL || inc2==NULL) return NULL;
	if (eval->collection == NULL) return NULL ;
	if (eval->result == NULL) return NULL ;
	if (eval->result->length != 1) return NULL ;

	params = (struct mask_params_s *)(eval->params) ;
	result = eval->result ;
	sc = eval->collection ;

	/* ensure we have two inputs which are evaluatable*/
	if (!sc_hasTwoInputs(sc)) return NULL  ;
	if (sc->input1->evaluator == NULL) return NULL ;
	if (sc->input2->evaluator == NULL) return NULL  ;

	/* ensure all the results have the same number of values */
	if ( !( (result->length != sc->input1->evaluator->result->length) &&
			(result->length != sc->input2->evaluator->result->length)) )
	{
		return NULL ;
	}

	if (inc1(sc->input1->inclusion, point)) {
		result1 = sc_evaluate(sc->input1, point) ;
		have_value = (result1 != NULL) ;
		/* populate the result with the result from the first input */
		if (have_value) {
			for (i=0; i<result->length; i++) {
				result->data[i] = result1->data[i] ;
			}
		}
	}

	if (!have_value && inc2(sc->input2, point)) {
		result2 = sc_evaluate(sc->input2, point) ;
		have_value = (result2 != NULL) ;
		/* populate the result with the result from the first input */
		if (have_value) {
			for (i=0; i<result->length; i++) {
				result->data[i] = result1->data[i] ;
			}
		}
	}

	if (!have_value) return NULL ;

	return result ;
}

VALUE *
first_value_evaluator(EVALUATOR *eval, LWPOINT *point)
{
	INCLUDE_FN *inc1 ;
	INCLUDE_FN *inc2 ;

	/* check for naughty users */
	if (eval == NULL || point == NULL) return NULL  ;
	if (eval->collection == NULL) return NULL ;
	if (!sc_hasTwoInputs(eval->collection)) return NULL ;

	inc1 = eval->collection->input1->inclusion->includes ;
	inc2 = eval->collection->input2->inclusion->includes ;

	return first_value_util(eval, point, inc1, inc2) ;
}

/** @} */ /* end of firsval_evaluator_i documentation group */

/** @} */ /* end of evaluator_i documentation group */

/*
 * Implementations of the "SpatialCollection" interface.
 * ===========================================================================
 */

/**
 * Single geometry wrapper: a "collection" backed by a single LWGEOM object.
 * This collection uses geometry_includes for inclusion and mask_evaluator
 *
 * @param geom the geometry to wrap
 * @
 */
SPATIAL_COLLECTION *
sc_create_geometry_wrapper(LWGEOM *geom, double inside, double outside)

