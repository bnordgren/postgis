#include "liblwgeom.h"
#include "spatial_collection.h"
#include "proj_api.h"
#include "string.h"


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
	if (geom->srid != point->srid) return 0;

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
 * \defgroup proj_wrap_include An includes implementation which reprojects the input coordinates.
 *
 * This is designed to be used as the "includes" component of a reprojecting
 * collection wrapper.
 * @{
 */

struct proj_wrap_inc_s {
	INCLUDES *wrapped ;
	projPJ source ;
	projPJ dest ;
};

int
projection_includes(INCLUDES *inc, LWPOINT *point)
{
	LWGEOM *point2d_g ;
	LWPOINT *point2d ;
	struct proj_wrap_inc_s *params ;
	int result ;

	if (inc == NULL || inc->params == NULL || point == NULL) return 0 ;
	params = (struct proj_wrap_inc_s *)(inc->params) ;

	/* copy the point so the outside world doesn't see it change. */
	point2d = lwpoint_make2d(point->srid,
			     lwpoint_get_x(point), lwpoint_get_y(point));

	/* note that this is a cast, not a copy */
	point2d_g = lwpoint_as_lwgeom(point2d) ;

	/* project the point */
	lwgeom_transform(point2d_g, params->dest, params->source) ;

	/* now call the wrapped object with the projected coordinates */
	result = params->wrapped->includes(params->wrapped, point2d) ;
	lwpoint_free(point2d) ;

	return result ;
}

INCLUDES *
sc_create_projection_includes(INCLUDES *wrapped, projPJ source, projPJ dest)
{
	struct proj_wrap_inc_s *params ;

	params = (struct proj_wrap_inc_s *)lwalloc(sizeof(struct proj_wrap_inc_s));
	if (params == NULL) return NULL ;

	/* initialize the parameters object */
	params->wrapped = wrapped ;
	params->source = source ;
	params->dest = dest ;

	/* create the includes object */
	return inc_create((PARAMETERS *)params, projection_includes, NULL);
}

void
sc_destroy_projection_includes(INCLUDES *dead)
{
	if (dead == NULL) return ;
	lwfree(dead->params) ;
	inc_destroy(dead) ;
}

/** @} */ /* end of proj_wrap_includes documentation group */

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

	if (inc == NULL || point == NULL ) return -1 ;
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
	struct mask_params_s *params ;
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

	return eval_create((PARAMETERS *)params, mask_evaluator, idx_fn, 1) ;
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
 * \defgroup proj_wrap_eval An evaluator implementation which reprojects the query point
 *
 * @{
 */

struct proj_wrap_eval_s {
	EVALUATOR *wrapped ;
	projPJ source ;
	projPJ dest ;
};

VALUE *
projection_evaluator(EVALUATOR *eval, LWPOINT *point)
{
	LWGEOM *point2d_g ;
	LWPOINT *point2d ;
	struct proj_wrap_eval_s *params ;
	VALUE *result ;

	if (eval == NULL || eval->params == NULL || point == NULL) return 0 ;
	params = (struct proj_wrap_eval_s *)(eval->params) ;

	/* copy the point so the outside world doesn't see it change. */
	point2d = lwpoint_make2d(point->srid,
			     lwpoint_get_x(point), lwpoint_get_y(point));

	/* note that this is a cast, not a copy */
	point2d_g = lwpoint_as_lwgeom(point2d) ;

	/* project the point */
	lwgeom_transform(point2d_g, params->dest, params->source) ;

	/* now call the wrapped object with the projected coordinates */
	result = params->wrapped->evaluate(params->wrapped, point2d) ;
	lwpoint_free(point2d) ;

	return result ;
}

EVALUATOR *
sc_create_projection_eval(EVALUATOR *wrapped, projPJ source, projPJ dest)
{
	struct proj_wrap_eval_s *params ;

	params = (struct proj_wrap_eval_s *)lwalloc(sizeof(struct proj_wrap_eval_s));
	if (params == NULL) return NULL ;

	/* initialize the parameters object */
	params->wrapped = wrapped ;
	params->source = source ;
	params->dest = dest ;

	/* create the includes object */
	return eval_create((PARAMETERS *)params, projection_evaluator, NULL, 1);
}

void
sc_destroy_projection_eval(EVALUATOR *dead)
{
	if (dead == NULL) return ;
	lwfree(dead->params) ;
	eval_destroy(dead) ;
}


/** @} */ /* end of proj_wrap_eval documentation group */


/**
 * \defgroup firstval_evaluator_i Evaluator to return a value from one of the collection's two inputs.
 *
 * This evaluator implements a simple default rule for determining the
 * value of a two-input collection. If the first input includes the query
 * point, then the value from the first input collection is returned.
 * Otherwise, the value from the second input collection is returned.
 *
 * This Evaluator requires that the two inputs have the same number of
 * bands. (e.g., the same "length" in their returned VALUEs). It also
 * assumes that the coordinate system for both input collections is the
 * same as the coordinate system for this object. For the evaluate and
 * includes functions, this means that the functions expect points in
 * the same srid. For evaluateIndex and includesIndex, this means that
 * the same grid definition applies to both inputs and the result.
 * @{
 */


/**
 * Utility function which embodies the logic to select a return value
 * from two inputs given a query point.
 *
 * @param eval pointer to the evaluator utility.
 * @param point the point for which a VALUE is needed.
 * @param index a true value indicates that point is an index into the raster
 */
VALUE *
first_value_util(EVALUATOR *eval,
		LWPOINT *point,
		int index)
{
	VALUE *result, *result_in1, *result_in2 ;
	SPATIAL_COLLECTION *sc ;
	INCLUDE_FN inc1, inc2 ;
	EVALUATOR_FN eval1, eval2 ;
	int have_value ;
	int i ;

	/* check for naughty users */
	if (eval == NULL || point == NULL) return NULL;
	if (eval->collection == NULL) return NULL ;
	if (eval->result == NULL) return NULL ;

	result = eval->result ;
	sc = eval->collection ;

	/* ensure we have two inputs which are evaluatable*/
	if (!sc_hasTwoInputs(sc)) return NULL  ;
	if (sc->input1->evaluator == NULL) return NULL ;
	if (sc->input2->evaluator == NULL) return NULL  ;
	if (sc->input1->inclusion == NULL) return NULL ;
	if (sc->input2->inclusion == NULL) return NULL  ;

	/* get the evaluator and includes functions for the inputs */
	if (!index) {
		eval1 = sc->input1->evaluator->evaluate ;
		eval2 = sc->input2->evaluator->evaluate ;
		inc1  = sc->input1->inclusion->includes ;
		inc2  = sc->input2->inclusion->includes ;
	} else {
		eval1 = sc->input1->evaluator->evaluateIndex ;
		eval2 = sc->input2->evaluator->evaluateIndex ;
		inc1  = sc->input1->inclusion->includesIndex ;
		inc2  = sc->input2->inclusion->includesIndex ;
	}

	/* ensure that the function pointers are all set */
	if ( (eval1 == NULL) || (eval2 == NULL) ||
		 (inc1 == NULL)  || (inc2 == NULL) )
		return NULL ;

	/* ensure all the results have the same number of values */
	if ( !((result->length == sc->input1->evaluator->result->length) &&
		   (result->length == sc->input2->evaluator->result->length)))
	{
		return NULL ;
	}

	have_value = 0 ;
	if (inc1(sc->input1->inclusion, point)) {
		result_in1 = eval1(sc->input1->evaluator, point) ;
		have_value = (result_in1 != NULL) ;
		/* populate the result with the result from the first input */
		if (have_value) {
			for (i=0; i<result->length; i++) {
				result->data[i] = result_in1->data[i] ;
			}
		}
	}

	if (!have_value && inc2(sc->input2->inclusion, point)) {
		result_in2 = eval2(sc->input2->evaluator, point) ;
		have_value = (result_in2 != NULL) ;
		/* populate the result with the result from the first input */
		if (have_value) {
			for (i=0; i<result->length; i++) {
				result->data[i] = result_in2->data[i] ;
			}
		}
	}

	if (!have_value) return NULL ;

	return result ;
}

/**
 * Implementation of EVALUATOR_FN for the case where point
 * represents some real world coordinate. It will rely on
 * the includes() and evaluate() functions of the inputs to
 * determine the result.
 */
VALUE *
first_value_evaluator(EVALUATOR *eval, LWPOINT *point)
{
	return first_value_util(eval, point, 0) ;
}

/**
 * Implementation of EVALUATOR_FN for the case where point
 * represents grid coordinates. It will rely on
 * the includesIndex() and evaluateIndex() functions of the inputs to
 * determine the result.
 */
VALUE *
first_value_evaluatorIndex(EVALUATOR *eval, LWPOINT *point)
{
	return first_value_util(eval, point, 1) ;
}

/**
 * Creates an instance of the Evaluator interface which
 * will return a value from one of the two inputs.
 */
EVALUATOR *
sc_create_first_value_evaluator(SPATIAL_COLLECTION *first, SPATIAL_COLLECTION *second)
{
	if (first == NULL || second == NULL) return NULL ;
	if (first->evaluator == NULL || second->evaluator == NULL) return NULL ;
	if (first->evaluator->result == NULL || second->evaluator->result == NULL) {
		return NULL  ;
	}

	if (first->evaluator->result->length != second->evaluator->result->length) {
		return NULL ;
	}
	return eval_create(NULL, first_value_evaluator, first_value_evaluatorIndex,
			first->evaluator->result->length) ;
}

/**
 * Destroys objects created by sc_create_first_value_evaluator. Do
 * not pass any other objects to this function.
 */
void
sc_destroy_first_value_evaluator(EVALUATOR *dead)
{
	eval_destroy(dead) ;
}

/** @} */ /* end of firstval_evaluator_i documentation group */

/** @} */ /* end of evaluator_i documentation group */

/**
 * \defgroup spatial_collection_i Implementations of the "SpatialCollection" interface.
 *
 * @{
 */

/**
 * \defgroup geo_wrap_collection A "collection" backed by a single LWGEOM object.
 * This collection uses geometry_includes for inclusion and mask_evaluator
 * for evaluation.
 * @{
 */

/**
 * Structure allows the geometry wrapper to take ownership of the
 * wrapped geometry, freeing it when done.
 */
struct geo_wrap_s {
	LWGEOM *geom ;
	int     owned ;
};

/**
 * Constructor to create a geometry wrapper collection object.
 *
 * @param geom the geometry to wrap
 * @param owned flag indicates that the collection should take ownership
 *              of the wrapped geometry object (e.g., geometry will be
 *              freed when collection is destroyed).
 * @param inside the "value" associated with the geometry's interior
 * @param outside the "value" associated with the geometry's exterior
 */
SPATIAL_COLLECTION *
sc_create_geometry_wrapper(LWGEOM *geom, int owned,
		double inside, double outside)
{
	INCLUDES *inclusion ;
	EVALUATOR *evaluator ;
	struct geo_wrap_s *params ;


	inclusion = sc_create_geometry_includes(geom) ;
	evaluator = sc_create_mask_evaluator(inside, outside, 0) ;
	params = (struct geo_wrap_s *)lwalloc(sizeof(struct geo_wrap_s)) ;

	if ( inclusion == NULL || evaluator==NULL || params == NULL )
	{
		if (inclusion != NULL ) sc_destroy_geometry_includes(inclusion) ;
		if (evaluator != NULL ) sc_destroy_mask_evaluator(evaluator) ;
		if (params != NULL ) lwfree(params) ;
		return NULL ;
	}

	/* record the geometry */
	params->geom = geom ;
	params->owned = owned ;

	/* ensure the geometry has an extent */
	lwgeom_add_bbox(geom);

	return sc_create(SPATIAL_ONLY, geom->srid, geom->bbox,
			params, inclusion, evaluator) ;
}

/**
 * Destroys an instance of the geometry wrapper spatial collection.
 * Only pass items in which were created with the sc_create_geometry_wrapper
 * function.
 */
void
sc_destroy_geometry_wrapper(SPATIAL_COLLECTION *dead)
{
	struct geo_wrap_s *params ;
	if (dead == NULL) return ;

	sc_destroy_geometry_includes(dead->inclusion) ;
	sc_destroy_mask_evaluator(dead->evaluator) ;
	if (dead->params != NULL) {
		params = (struct geo_wrap_s *)(dead->params) ;
		if ((params->geom != NULL) && params->owned) {
			lwgeom_free(params->geom) ;
		}
		lwfree(params) ;
	}
	sc_destroy(dead) ;
}

/** @} */  /* end of geo_wrap_collection documentation group */

/**
 * \defgroup proj_wrap_collection Collection wrapper which reprojects coordinates.
 *
 * This collection implementation intercepts and reprojects the
 * real-world coordinates passed to the evaluate and includes methods.
 * At creation time, the user declares the projection in which they
 * will supply coordinates. Whenever the evaluate or includes methods
 * are called, the point is projected from the declared projection
 * to the projection expected by the collection which this object is
 * wrapping.
 *
 * @{
 */

struct proj_wrap_s {
	SPATIAL_COLLECTION *wrapped ;
	projPJ source ;
	projPJ dest ;
};

LWPOLY *
gbox_to_lwpoly(GBOX *bbox)
{
	POINTARRAY **rings ;
	POINTARRAY *pts ;
	POINT4D     p4d ;

    rings = (POINTARRAY **) lwalloc(sizeof (POINTARRAY*));
    if (!rings) {
        return NULL ;
    }
    rings[0] = ptarray_construct(0, 0, 5);
    /* TODO: handle error on ptarray construction */
    /* XXX jorgearevalo: the error conditions aren't managed in ptarray_construct */
    if (!rings[0]) {
        return NULL;
    }
    pts = rings[0];

    /* first corner */
    p4d.x = bbox->xmin ;
    p4d.y = bbox->ymin ;
    ptarray_set_point4d(pts, 0, &p4d);
    ptarray_set_point4d(pts, 4, &p4d); /* needed for closing it? */

    /* second corner */
    p4d.x = bbox->xmin ;
    p4d.y = bbox->ymax ;
    ptarray_set_point4d(pts, 1, &p4d);

    /* third corner */
    p4d.x = bbox->xmax ;
    p4d.y = bbox->ymax ;
    ptarray_set_point4d(pts, 2, &p4d);

    /* fourth corner */
    p4d.x = bbox->xmax ;
    p4d.y = bbox->ymin ;
    ptarray_set_point4d(pts, 3, &p4d);

    /* make the polygon */
    return lwpoly_construct(SRID_UNKNOWN, 0, 1, rings);
}

void
project_gbox(GBOX *bbox, projPJ from, projPJ to)
{
	LWPOLY *extent ;

	/* make a polygon out of the gbox */
	extent = gbox_to_lwpoly(bbox) ;
	if (extent == NULL) return ;

    /* project it and recalculate the mins and maxes */
    lwgeom_transform(lwpoly_as_lwgeom(extent), from, to) ;
    lwgeom_drop_bbox(lwpoly_as_lwgeom(extent)) ;
    lwgeom_add_bbox(lwpoly_as_lwgeom(extent)) ;

    /* copy this over the input bbox */
	memcpy(bbox, extent->bbox, sizeof(GBOX)) ;

	lwpoly_free(extent) ;
}

SPATIAL_COLLECTION *
sc_create_projection_wrapper(SPATIAL_COLLECTION *wrapped,
		                     int32_t desired_srid,
		                     projPJ wrapped_proj, projPJ desired_proj )
{
	INCLUDES *inclusion ;
	EVALUATOR *evaluator ;
	GBOX extent ;
	struct proj_wrap_s *params ;

	if (wrapped == NULL) return NULL ;

	params = (struct proj_wrap_s *)lwalloc(sizeof(struct proj_wrap_s)) ;
	if (params == NULL) return NULL ;

	/* initialize the parameters */
	params->wrapped = wrapped ;
	params->source  = desired_proj ;
	params->dest    = wrapped_proj ;

	inclusion = sc_create_projection_includes(wrapped->inclusion,
			         desired_proj, wrapped_proj) ;
	evaluator = sc_create_projection_eval(wrapped->evaluator,
			         desired_proj, wrapped_proj) ;

	if ( inclusion == NULL || evaluator==NULL )
	{
		if (inclusion != NULL ) sc_destroy_projection_includes(inclusion) ;
		if (evaluator != NULL ) sc_destroy_projection_eval(evaluator) ;
		return NULL ;
	}

	/* calculate the projected extent */
	memcpy(&extent, &(wrapped->extent), sizeof(GBOX)) ;
	project_gbox(&extent, wrapped_proj, desired_proj) ;

	return sc_create(wrapped->type, desired_srid, &extent,
			         (PARAMETERS *)params, inclusion, evaluator) ;
}

void
sc_destroy_projection_wrapper(SPATIAL_COLLECTION *dead)
{
	if (dead == NULL) return ;

	if (dead->params != NULL) {
		lwfree(dead->params) ;
	}
	sc_destroy_projection_includes(dead->inclusion) ;
	sc_destroy_projection_eval(dead->evaluator) ;
	sc_destroy(dead) ;
}



/** @} */  /* end of proj_wrap_collection documentation group */

/**
 * \defgroup relation_collection Creates a collection formed by the spatial relationship of two input collections.
 * @{
 */

/**
 * Calculates the extent of two collections related by the intersection
 * operation.
 */
GBOX *relation_env_intersection(LWPOLY *r1, LWPOLY *r2)
{
	LWGEOM *env ;
	GBOX *bounds ;


	env = lwgeom_intersection(lwpoly_as_lwgeom(r1),
			                   lwpoly_as_lwgeom(r2)) ;
	lwgeom_add_bbox(env) ;
	bounds = gbox_copy(env->bbox);
	lwgeom_free(env) ;
	return bounds ;
}

/**
 * Calculates the extent of two collections related by the difference
 * operation.
 */
GBOX *relation_env_difference(LWPOLY *r1, LWPOLY *r2)
{
	LWGEOM *env ;
	GBOX *bounds ;

	env = lwgeom_difference(lwpoly_as_lwgeom(r1),
			                   lwpoly_as_lwgeom(r2)) ;
	lwgeom_add_bbox(env) ;
	bounds = gbox_copy(env->bbox);
	lwgeom_free(env) ;
	return bounds ;
}

/**
 * Calculates the extent of two collections related by the union
 * operation.
 */
GBOX *relation_env_union(LWPOLY *r1, LWPOLY *r2)
{
	LWGEOM *env ;
	GBOX *bounds ;

	env = lwgeom_union(lwpoly_as_lwgeom(r1),
			                   lwpoly_as_lwgeom(r2)) ;
	lwgeom_add_bbox(env) ;
	bounds = gbox_copy(env->bbox);
	lwgeom_free(env) ;
	return bounds ;
}

/**
 * Calculates the extent of two collections related by the symmetric
 * difference operation.
 */
GBOX *relation_env_symdifference(LWPOLY *r1, LWPOLY *r2)
{
	LWGEOM *env ;
	GBOX *bounds ;

	env = lwgeom_symdifference(lwpoly_as_lwgeom(r1),
			                   lwpoly_as_lwgeom(r2)) ;
	lwgeom_add_bbox(env) ;
	bounds = gbox_copy(env->bbox);
	lwgeom_free(env) ;
	return bounds ;
}

/**
 * \defgroup relation_op Relates two input collections in the same projection.
 *
 * @{
 */

/**
 * performs the heavy lifting for #sc_create_relation_op and
 * #sc_create_relation_op_proj.
 */
static SPATIAL_COLLECTION *
sc_create_relation_op_internal(COLLECTION_TYPE t,
		              SPATIAL_COLLECTION *sc1,
		              SPATIAL_COLLECTION *sc2,
		              PARAMETERS *params,
		              ENVELOPE_PREP_OP env_fn,
		              RELATION_FN inc_fn,
		              EVALUATOR *eval)
{
	SPATIAL_COLLECTION *result ;
	INCLUDES *inc ;
	LWPOLY *sc1_outline ;
	LWPOLY *sc2_outline ;
	GBOX *result_env ;

	/* sanity checks on inputs */
	if (sc1 == NULL || sc2==NULL) return NULL ;
	if (env_fn == NULL || inc_fn == NULL) return NULL ;
	if (eval == NULL) return NULL ;

	/* something which should print out an error */
	if (sc1->srid != sc2->srid) return NULL ;

	/* create the INCLUDES given the desired function */
	inc = sc_create_relation_includes(inc_fn) ;
	if (inc == NULL) return NULL ;

	/* calculate an approximate envelope for the result */
	sc1_outline = gbox_to_lwpoly(&(sc1->extent)) ;
	sc2_outline = gbox_to_lwpoly(&(sc2->extent)) ;
	if (sc1_outline == NULL || sc2_outline == NULL) {
		sc_destroy_relation_includes(inc) ;
		if (sc1_outline != NULL) lwpoly_free(sc1_outline) ;
		if (sc2_outline != NULL) lwpoly_free(sc2_outline) ;
		return NULL ;
	}
	result_env = env_fn(sc1_outline, sc2_outline) ;
	lwpoly_free(sc1_outline) ;
	lwpoly_free(sc2_outline) ;
	if (result_env == NULL) {
		sc_destroy_relation_includes(inc) ;
		return NULL ;
	}

	/* create the two-input collection */
	result = sc_twoinput_create(t, params, result_env, inc, eval, sc1, sc2) ;
	if (result == NULL) {
		sc_destroy_relation_includes(inc) ;
		return NULL ;
	}

	return result ;
}

/**
 * A collection implementation encapsulating the spatial relationship
 * of two input collections. The user may specify an existing
 * #ENVELOPE_PREP_OP to calculate the spatial extent of the result.
 * A #RELATION_FN may also be provided to determine whether an
 * individual point is included in the result. An arbitrary
 * #EVALUATOR may also be provided (which must be compatible with
 * the inputs).
 *
 * The result and the two inputs must be in the same projection
 * (same srid).
 *
 * This constructor is probably not what you want, unless you have
 * written a custom #ENVELOPE_PREP_OP or #RELATION_FN (or need to
 * specify some nonstandard combination.) See the #sc_create_sync_relation_op
 * for a constructor which synchronizes the #ENVELOPE_PREP_OP and
 * #RELATION_FN to an enumerated #RELATION_TYPE.
 *
 * When the result is expressed in a different projection than
 * one or both inputs, see #sc_create_relation_op_proj or
 * #sc_create_sync_relation_op_proj.
 *
 * @returns the resultant collection
 * @param t type of the collection (#SPATIAL_ONLY or #SPATIAL_PLUS_VALUE).
 * @param sc1 the first input collection
 * @param sc2 the second input collection
 * @param env_fn an #ENVELOPE_PREP_OP to calculate the extent of the result
 *               given sc1 and sc2
 * @param inc_fn a #RELATION_FN to calculate whether an individual point
 *               is included in the result, given sc1 and sc2
 * @param eval   an #EVALUATOR which provides a resultant value given
 *               the values of sc1 and/or sc2. May be NULL if t is
 *               #SPATIAL_ONLY
 */
SPATIAL_COLLECTION *
sc_create_relation_op(COLLECTION_TYPE t,
		              SPATIAL_COLLECTION *sc1,
		              SPATIAL_COLLECTION *sc2,
		              ENVELOPE_PREP_OP env_fn,
		              RELATION_FN inc_fn,
		              EVALUATOR *eval)
{
	return sc_create_relation_op_internal(t, sc1, sc2, NULL,
			env_fn, inc_fn, eval) ;
}



/**
 * A collection implementation encapsulating the spatial relationship
 * of two input collections. The user specifies the relation by supplying
 * a #RELATION_TYPE value.  An arbitrary
 * #EVALUATOR may also be provided (which must be compatible with
 * the inputs).
 *
 * The result and the two inputs must be in the same projection
 * (same srid).
 *
 * When the result is expressed in a different projection than
 * one or both inputs, see #sc_create_relation_op_proj or
 * #sc_create_sync_relation_op_proj.
 *
 * @returns the resultant collection
 * @param t type of the collection (#SPATIAL_ONLY or #SPATIAL_PLUS_VALUE).
 * @param sc1 the first input collection
 * @param sc2 the second input collection
 * @param relation specifies a supported, predefined spatial relationship
 * @param eval   an #EVALUATOR which provides a resultant value given
 *               the values of sc1 and/or sc2. May be NULL if t is
 *               #SPATIAL_ONLY
 */
SPATIAL_COLLECTION *
sc_create_sync_relation_op(COLLECTION_TYPE t,
		              SPATIAL_COLLECTION *sc1,
		              SPATIAL_COLLECTION *sc2,
		              RELATION_TYPE relation,
		              EVALUATOR *eval)
{
	RELATION_FN      inc_fn ;
	ENVELOPE_PREP_OP env_fn ;

	if (sc1 == NULL || sc2 == NULL) return NULL ;

	/*
	 * get the envelope and includes functions
	 * (and make sure they correspond).
	 */
	inc_fn = sc_get_relation_fn(relation) ;
	env_fn = sc_get_envelope_fn(relation) ;
	if (inc_fn == NULL || env_fn == NULL) return NULL ;

	return sc_create_relation_op(t, sc1, sc2,
			env_fn, inc_fn, eval) ;
}


void
sc_destroy_relation_op(SPATIAL_COLLECTION *dead)
{
	if (dead != NULL) {
		if (dead->inclusion != NULL) {
			sc_destroy_relation_includes(dead->inclusion) ;
		}
		sc_twoinput_destroy(dead) ;
	}
}

/** @} */ /* end of the relation_op documentation group */

/**
 * \defgroup relation_proj Wraps each input with a projection wrapper, and installs the specified relationship functions.
 *
 * @{
 */

struct relation_op_proj_s {
	SPATIAL_COLLECTION *sc1_wrap ;
	SPATIAL_COLLECTION *sc2_wrap ;
};


/**
 * A collection implementation encapsulating the spatial relationship
 * of two input collections. The user may specify an existing
 * #ENVELOPE_PREP_OP to calculate the spatial extent of the result.
 * A #RELATION_FN may also be provided to determine whether an
 * individual point is included in the result. An arbitrary
 * #EVALUATOR may also be provided (which must be compatible with
 * the inputs).
 *
 * Projection information must be specified for the result and
 * both inputs. If necessary, each input is wrapped with a
 * projection wrapper into the destination projection.
 *
 * This constructor is probably not what you want, unless you have
 * written a custom #ENVELOPE_PREP_OP or #RELATION_FN (or need to
 * specify some nonstandard combination.) See the #sc_create_sync_relation_op
 * for a constructor which synchronizes the #ENVELOPE_PREP_OP and
 * #RELATION_FN to an enumerated #RELATION_TYPE.
 *
 *
 * When the result and the two inputs are known to be in the same
 * projection, see #sc_create_relation_op or
 * #sc_create_sync_relation_op.
 *
 * @returns the resultant collection
 * @param t type of the collection (#SPATIAL_ONLY or #SPATIAL_PLUS_VALUE).
 * @param sc1 the first input collection
 * @param sc2 the second input collection
 * @param proj_sc1 projection information for first input
 * @param proj_sc2 projection information for second input
 * @param srid     srid of the result
 * @param proj_dest projection information for the result
 * @param env_fn an #ENVELOPE_PREP_OP to calculate the extent of the result
 *               given sc1 and sc2
 * @param inc_fn a #RELATION_FN to calculate whether an individual point
 *               is included in the result, given sc1 and sc2
 * @param eval   an #EVALUATOR which provides a resultant value given
 *               the values of sc1 and/or sc2. May be NULL if t is
 *               #SPATIAL_ONLY
 */

SPATIAL_COLLECTION *
sc_create_relation_op_proj(COLLECTION_TYPE t,
		              SPATIAL_COLLECTION *sc1,
		              SPATIAL_COLLECTION *sc2,
		              projPJ proj_sc1, projPJ proj_sc2,
		              int32_t srid, projPJ proj_dest,
		              ENVELOPE_PREP_OP env_fn,
		              RELATION_FN inc_fn,
		              EVALUATOR *eval)
{
	SPATIAL_COLLECTION *sc1_used ;
	SPATIAL_COLLECTION *sc2_used ;
	struct relation_op_proj_s *params ;

	if (sc1 == NULL || sc2 == NULL) return NULL ;
	if (proj_sc1 == NULL || proj_sc2 == NULL || proj_dest == NULL) return NULL ;

	params = (struct relation_op_proj_s *)lwalloc(sizeof(struct relation_op_proj_s));
	if (params==NULL) return NULL ;
	params->sc1_wrap = NULL ;
	params->sc2_wrap = NULL ;

	/* wrap sc1 if necessary */
	sc1_used = sc1 ;
	if (sc1->srid != srid) {
		params->sc1_wrap =
				sc_create_projection_wrapper(sc1, srid, proj_sc1, proj_dest) ;
		if (params->sc1_wrap == NULL) {
			lwfree(params) ;
			return NULL ;
		}
		sc1_used = params->sc1_wrap ;
	}

	/* wrap sc2 if necessary */
	sc2_used = sc2 ;
	if (sc2->srid != srid) {
		params->sc2_wrap =
				sc_create_projection_wrapper(sc2, srid, proj_sc2, proj_dest) ;
		if (params->sc2_wrap == NULL) {
			if (params->sc1_wrap != NULL) {
				sc_destroy_projection_wrapper(params->sc1_wrap) ;
			}
			lwfree(params) ;
			return NULL ;
		}
		sc2_used = params->sc2_wrap ;
	}

	return sc_create_relation_op_internal(t, sc1_used, sc2_used,
			params, env_fn, inc_fn, eval) ;
}

/**
 * A collection implementation encapsulating the spatial relationship
 * of two input collections. The user may specifies a supported,
 * predefined spatial operation using the #RELATION_TYPE enumeration.
 * An arbitrary
 * #EVALUATOR may also be provided (which must be compatible with
 * the inputs).
 *
 * Projection information must be specified for
 * both inputs. The projection of the first input is used as the
 * projection of the result.
 *
 * When the result and the two inputs are known to be in the same
 * projection, see #sc_create_relation_op or
 * #sc_create_sync_relation_op.
 *
 * @returns the resultant collection
 * @param t type of the collection (#SPATIAL_ONLY or #SPATIAL_PLUS_VALUE).
 * @param sc1 the first input collection
 * @param sc2 the second input collection
 * @param proj_sc1 projection information for first input
 * @param proj_sc2 projection information for second input
 * @param relation the desired predefined spatial relationship
 * @param eval   an #EVALUATOR which provides a resultant value given
 *               the values of sc1 and/or sc2. May be NULL if t is
 *               #SPATIAL_ONLY
 */
SPATIAL_COLLECTION *
sc_create_sync_relation_op_proj(COLLECTION_TYPE t,
		              SPATIAL_COLLECTION *sc1,
		              SPATIAL_COLLECTION *sc2,
		              projPJ proj_sc1, projPJ proj_sc2,
		              int dest_srid, projPJ proj_dest,
		              RELATION_TYPE relation,
		              EVALUATOR *eval)
{
	RELATION_FN      inc_fn ;
	ENVELOPE_PREP_OP env_fn ;

	if (sc1 == NULL || sc2 == NULL) return NULL ;

	/*
	 * get the envelope and includes functions
	 * (and make sure they correspond).
	 */
	inc_fn = sc_get_relation_fn(relation) ;
	env_fn = sc_get_envelope_fn(relation) ;
	if (inc_fn == NULL || env_fn == NULL) return NULL ;

	return sc_create_relation_op_proj(t, sc1, sc2,
			proj_sc1, proj_sc2,
			dest_srid, proj_dest,
			env_fn, inc_fn, eval) ;
}

void
sc_destroy_relation_op_proj(SPATIAL_COLLECTION *dead)
{
	if (dead != NULL) {
		if (dead->params != NULL) {
			struct relation_op_proj_s *p ;
			p = (struct relation_op_proj_s *)(dead->params) ;
			if (p->sc1_wrap != NULL) {
				sc_destroy_projection_wrapper(p->sc1_wrap) ;
			}
			if (p->sc2_wrap != NULL) {
				sc_destroy_projection_wrapper(p->sc2_wrap) ;
			}
		}
	}
	sc_destroy_relation_op(dead) ;
}

/** @} */ /* end of the relation_proj documentation group */


/** @} */ /* end of relation_collection documentation group */

/** @} */  /* end of spatial_collection_i documentation group */


/**
 * Converts a #RELATION_TYPE value to a #RELATION_FN.
 * Returns NULL if the value is not recognized.
 */
RELATION_FN
sc_get_relation_fn(RELATION_TYPE relation)
{
	RELATION_FN result ;

	result = NULL ;
	switch (relation) {
	case INTERSECTION :
		result = relation_intersection ;
		break ;
	case UNION :
		result = relation_union ;
		break ;
	case DIFFERENCE :
		result = relation_difference ;
		break ;
	case SYMDIFFERENCE :
		result = relation_symdifference ;
		break ;
	default:
		; /* leave it null */
	}

	return result ;
}

/**
 * Converts a #RELATION_TYPE value to an
 * #ENVELOPE_PREP_OP function pointer. Returns
 * NULL if the value is not recognized.
 */
ENVELOPE_PREP_OP
sc_get_envelope_fn(RELATION_TYPE relation)
{
	ENVELOPE_PREP_OP result ;

	result = NULL ;
	switch (relation) {
	case INTERSECTION :
		result = relation_env_intersection ;
		break ;
	case UNION :
		result = relation_env_union ;
		break ;
	case DIFFERENCE :
		result = relation_env_difference ;
		break ;
	case SYMDIFFERENCE :
		result = relation_env_symdifference ;
		break ;
	default:
		; /* leave it null */
	}

	return result ;
}