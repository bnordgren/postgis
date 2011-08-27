#include "sc_raster.h"
#include "rt_api.h"
#include "gdal.h"
#include "math.h"


/**
 * \addtogroup includes_i Implementations of the "includes" interface.
 * @{
 */

/**
 * \defgroup rt_env_inc Considers all grid cells to be included.
 *
 * @{
 * This implementation of the includes interface will return true if
 * the provided location is somewhere on the grid defined by the
 * raster. The Nodata value, if it exists, is ignored. This implementation
 * supports both the includesIndex() and includes() methods.
 */


struct raster_inc_s {
	rt_raster source ;
	int band_num ;
};


/**
 * Considers an index to be included iff the point is between
 * (0,0) and (width-1, height-1).
 * @param inc pointer to the includes structure
 * @param point the point to evaluate
 * @return true if point is within the bounds of the raster.
 */
int
raster_envelope_includes_idx(INCLUDES *inc, LWPOINT *point)
{
	struct raster_inc_s *p ;
	uint16_t width ;
	uint16_t height ;
	POINT2D  indices ;

	if (inc == NULL || point == NULL) return 0 ;
	if (inc->params == NULL) return 0 ;

	p = (struct raster_inc_s *)(inc->params);
	width = rt_raster_get_width(p->source) ;
	height = rt_raster_get_height(p->source) ;

	lwpoint_getPoint2d_p(point, &indices) ;
	return (indices.x >= 0) && (indices.x < width) &&
		   (indices.y >= 0) && (indices.y < height) ;
}

/**
 * Considers a point to be included iff it is expressed in the
 * same coordinate system as the raster. The point is then transformed
 * into image indices and checked using the includesIndex() function. This
 * is a means of wrapping the index function.
 * @param inc pointer to the includes structure
 * @param point the point to evaluate (must be in same srid as the raster)
 * @return true if point is within the bounds of the raster.
 */
int
raster_idxwrap_includes(INCLUDES *inc, LWPOINT *point)
{
	struct raster_inc_s *p ;
	POINT2D world ;
	POINT2D indices ;
	LWPOINT *idx_point ;
	int included ;

	if (inc == NULL || point == NULL) return 0 ;
	if (inc->params == NULL) return 0 ;
	if (inc->includesIndex == NULL) {
		rterror("raster_idxwrap_includes: the includesIndex function must be specified.") ;
		return 0 ;
	}

	/*
	 * check that the evaluation point and the raster are expressed
	 * in the same coordinate system.
	 */
	p = (struct raster_inc_s *)(inc->params);
	if (p->source == NULL ) return 0 ;
	if (rt_raster_get_srid(p->source) != point->srid) {
		rterror("rt_idxwrap_includes: the point must have the same srid as the raster");
		return 0 ;
	}

	/* convert world coordinates to raster indices */
	lwpoint_getPoint2d_p(point, &world) ;
	rt_raster_geopoint_to_cell(p->source,
			              world.x, world.y,
			              &(indices.x),&(indices.y)) ;

	/* check whether the index is "in range" */
	idx_point = lwpoint_make2d(SRID_UNKNOWN, indices.x, indices.y) ;
	included = inc->includesIndex(inc, idx_point) ;
	lwpoint_free(idx_point) ;

	return included;
}

/**
 * Creates an instance of a simple and fast implementation of the
 * Includes interface. This implementation checks only to see that the
 * supplied coordinates are inside the grid specified by the raster. It
 * does not check the raster for nodata values. This is appropriate for
 * cases where the raster does not have a nodata value or when only the
 * bounds are interesting.
 */
INCLUDES *
sc_create_raster_env_includes(rt_raster raster)
{
	struct raster_inc_s *params ;
	INCLUDES *inc ;

	params = (struct raster_inc_s *)(rtalloc(sizeof(struct raster_inc_s)));
	if (params == NULL) return NULL ;

	params->source = raster ;
	inc = inc_create(params, raster_idxwrap_includes, raster_envelope_includes_idx);
	if (inc == NULL) {
		rtdealloc(params) ;
	}

	return inc ;
}
/**
 * Destroys an instance of the Includes interface created with
 * sc_create_raster_env_includes(). Do not pass any Includes
 * object created by another function.
 */
void
sc_destroy_raster_env_includes(INCLUDES *dead)
{
	if (dead != NULL) {
		if (dead->params != NULL) {
			rtdealloc(dead->params) ;
			inc_destroy(dead) ;
		}
	}
}

/** @} */ /* end of rt_env_inc documentation group */


/**
 * \defgroup rt_nodata_inc An includes implementation which accounts for "nodata".
 *
 * @{
 * This implementation of the includes interface only considers a
 * point to be included if it occurs within the envelope of the
 * raster and if the value of the associated pixel in the specified band is
 * something other than "nodata". Nearest neighbor resampling is
 * performed on the evaluation point to select a pixel.
 */

/**
 * After rounding the supplied coordinates (assumed to be indices) to the
 * nearest integer, the new coordinates are checked to ensure that they are
 * within the raster. If so, the raster (not the collection) is evaluated at
 * that grid cell and the value is compared to the nodata value.
 *
 * @param inc the includes structure
 * @param point the point to evaluate (in the image indices coordinate system).
 * @return true if the raster contains data at the pixel nearest the specified point.
 */
int
raster_nodata_includes_idx(INCLUDES *inc, LWPOINT *point)
{
	struct raster_inc_s *p ;
	rt_band band ;
	POINT2D rounded ;
	LWPOINT *rounded_lwp ;
	int included ;

	if (inc == NULL || point == NULL) return 0  ;
	if (inc->params == NULL) return 0  ;

	p = (struct raster_inc_s *)(inc->params) ;
	if (p->source == NULL) return 0  ;

	/* round the coordinates */
	lwpoint_getPoint2d_p(point, &rounded) ;
	rounded.x = round(rounded.x) ;
	rounded.y = round(rounded.y) ;
	rounded_lwp = lwpoint_make2d(SRID_UNKNOWN, rounded.x, rounded.y) ;

	/* check the rounded point is in the envelope */
	included = raster_envelope_includes_idx(inc, rounded_lwp) ;

	/* now check the value at the point, if necessary */
	if (included) {

		/* get the band to test for nodata */
		band = rt_raster_get_band(p->source, p->band_num) ;
		if (band == NULL) {
			rterror("raster_nodata_includes_idx: raster does not contain band %d\n",
					p->band_num) ;
			lwpoint_free(rounded_lwp) ;
			return 0 ;
		}

		if (rt_band_get_hasnodata_flag(band)) {
			double raster_val ;

			/* get the value at the pixel. */
			rt_band_get_pixel(band, (uint16_t)(rounded.x), (uint16_t)(rounded.y), &raster_val);

			/* compare with nodata */
			included = (raster_val != rt_band_get_nodata(band)) ;
		}
	}

	lwpoint_free(rounded_lwp) ;

	return included;
}

/**
 * Creates an instance of the Includes interface which is backed by a raster
 * having a band containing nodata values. The includes() and includesIndex()
 * functions will return true iff the specified point is within the envelope
 * of the raster and the value of the nearest pixel is something other than
 * nodata.
 * @param r the raster which will "back" this instance of the Includes interface
 * @param band the band to examine for nodata values.
 * @return a new instance of the Includes interface which examines a raster for
 *        nodata values.
 */
INCLUDES *
sc_create_raster_nodata_includes(rt_raster r, int band)
{
	INCLUDES *inc ;
	struct raster_inc_s *params ;
	rt_band test_band ;

	if (r == NULL) return NULL ;

	/* check that the raster has the specified band */
	test_band = rt_raster_get_band(r, band) ;
	if (test_band == NULL) {
		rterror("sc_create_raster_nodata_includes: invalid band %d\n",band) ;
		return NULL ;
	}

	params = (struct raster_inc_s *)rtalloc(sizeof (struct raster_inc_s));
	if (params == NULL) return NULL ;
	params->source = r ;
	params->band_num = band ;

	inc = inc_create(params, raster_idxwrap_includes, raster_nodata_includes_idx) ;

	return inc ;
}

/**
 * Destroys an instance of the Includes interface created with
 * sc_create_raster_nodata_includes(). Do not pass any Includes
 * object created by another function.
 */
void
sc_destroy_raster_nodata_includes(INCLUDES *dead)
{
	if (dead != NULL) {
		if (dead->params != NULL) {
			rtdealloc(dead->params) ;
			inc_destroy(dead) ;
		}
	}
}

/** @} */ /* end of documentation group rt_nodata_inc */

/** @} */ /* end of documentation group includes_i */

/**
 * \addtogroup evaluator_i Implementations of the "evaluator" interface.
 * @{
 */

/**
 * \defgroup raster_eval Returns the value(s) of the specified band(s).
 *
 * @{
 * This evaluator allows the user to specify one or more raster bands
 * to use when generating VALUEs. The default is to use all of the bands
 * in order.
 */

struct raster_eval_s {
	rt_raster source ;
	int num_bands ;
	int *bands ;
};

VALUE *
raster_bands_evaluate_idx(EVALUATOR *eval, LWPOINT *point)
{
	POINT2D indices ;
	struct raster_eval_s *p ;
	int goodPt;
	int i ;

	if (eval == NULL || point == NULL) return NULL ;
	if (eval->params == NULL) return NULL ;
	if (eval->result == NULL) return NULL ;

	p = (struct raster_eval_s *)(eval->params) ;
	if (p->source == NULL || p->bands == NULL) return NULL ;
	if (p->num_bands != eval->result->length) return NULL ;

	/* check bounds against raster */
	lwpoint_getPoint2d_p(point, &indices) ;
	goodPt = (indices.x >= 0) && (indices.x < rt_raster_get_width(p->source)) &&
			 (indices.y >= 0) && (indices.y < rt_raster_get_height(p->source)) ;
	if (!goodPt) return NULL ;

	/* loop over all the specified bands and get the same pixel in each one */
	for (i=0; i<p->num_bands; i++) {
		rt_band band ;
		double band_val ;

		band = rt_raster_get_band(p->source, p->bands[i]) ;
		if (band == NULL) {
			rterror("raster_bands_evaluate_idx: cannot get band %d", p->bands[i]) ;
			return NULL ;
		}

		rt_band_get_pixel(band, indices.x, indices.y, &band_val) ;

		eval->result->data[i] = band_val ;
	}

	return eval->result ;
}

VALUE *
raster_idxwrap_evaluate(EVALUATOR *eval, LWPOINT *point)
{
	struct raster_eval_s *p ;
	POINT2D world ;
	POINT2D indices ;
	LWPOINT *idx_point ;
	VALUE *result ;

	if (eval == NULL || point == NULL) return 0 ;
	if (eval->params == NULL) return 0 ;
	if (eval->evaluateIndex == NULL) {
		rterror("raster_idxwrap_evaluate: the evaluateIndex function must be specified.") ;
		return 0 ;
	}

	/*
	 * check that the evaluation point and the raster are expressed
	 * in the same coordinate system.
	 */
	p = (struct raster_eval_s *)(eval->params);
	if (p->source == NULL ) return 0 ;
	if (rt_raster_get_srid(p->source) != point->srid) {
		rterror("rt_idxwrap_evaluate: the point must have the same srid as the raster");
		return 0 ;
	}

	/* convert world coordinates to raster indices */
	lwpoint_getPoint2d_p(point, &world) ;
	rt_raster_geopoint_to_cell(p->source,
			              world.x, world.y,
			              &(indices.x),&(indices.y)) ;

	/* evaluate the given index */
	idx_point = lwpoint_make2d(SRID_UNKNOWN, indices.x, indices.y) ;
	result = eval->evaluateIndex(eval, idx_point) ;
	lwpoint_free(idx_point) ;

	return result;
}

EVALUATOR *
sc_create_raster_bands_evaluator(rt_raster raster, int *bands, int num_bands)
{
	rt_band test_band ;
	struct raster_eval_s *params ;
	EVALUATOR *eval ;
	int i ;

	if (raster == NULL) return NULL ;

	if (bands == NULL) {
		/* the default is to pull all the bands in order */
		num_bands = rt_raster_get_num_bands(raster) ;
	} else {
		/* otherwise we try to do what they ask */
		for (i=0; i<num_bands; i++) {
			test_band = rt_raster_get_band(raster, bands[i]) ;
			if (test_band == NULL) {
				rterror("sc_create_raster_bands_evaluator: cannot load band %d\n",bands[i]);
				return NULL ;
			}
		}
	}

	/** allocate space for the parameters. */
	params = (struct raster_eval_s *)rtalloc(sizeof(struct raster_eval_s)) ;
	if (params == NULL) return NULL ;
	params->bands = (int*)rtalloc(sizeof(int)*num_bands) ;
	if (params->bands == NULL) {
		rtdealloc(params) ;
		return NULL ;
	}

	/* initialize the parameter structure */
	params->num_bands = num_bands ;
	params->source = raster ;
	for (i=0; i<num_bands; i++) {
		params->bands[i] = (bands==NULL) ? i : bands[i] ;
	}

	eval = eval_create(params,
			raster_idxwrap_evaluate, raster_bands_evaluate_idx,
			num_bands) ;

	if (eval == NULL) {
		rtdealloc(params->bands) ;
		rtdealloc(params) ;
		return NULL ;
	}

	return eval ;
}

void
sc_destroy_raster_bands_evaluator(EVALUATOR *dead)
{
	if (dead != NULL) {
		if (dead->params != NULL) {
			struct raster_eval_s *dead_params ;

			dead_params = (struct raster_eval_s *)(dead->params) ;
			if (dead_params->bands != NULL) {
				rtdealloc(dead_params->bands) ;
			}
			rtdealloc(dead_params) ;
		}
		eval_destroy(dead) ;
	}
}

/** @} */ /* end of documentation group raster_eval */

/** @} */ /* end of documentation group evaluator_i */

/**
 * \addtogroup spatial_collection_i Implementations of the "SpatialCollection" interface.
 * @{
 */

SPATIAL_COLLECTION *
sc_create_raster_wrapper(rt_raster raster, int *bands, int num_bands)
{
	INCLUDES *inc ;
	EVALUATOR *eval ;
	LWPOLY *outline ;
	SPATIAL_COLLECTION *sc ;

	if (raster==NULL) return NULL ;

	inc = sc_create_raster_env_includes(raster) ;
	if (inc == NULL) return NULL ;

	eval = sc_create_raster_bands_evaluator(raster, bands, num_bands) ;
	if (eval == NULL) {
		sc_destroy_raster_env_includes(inc) ;
		return NULL ;
	}

	outline = rt_raster_get_convex_hull(raster) ;
	lwgeom_add_bbox(lwpoly_as_lwgeom(outline)) ;

	sc = sc_create(SPATIAL_PLUS_VALUE,
			rt_raster_get_srid(raster),
			outline->bbox,
			NULL, inc, eval) ;

	if (sc == NULL) {
		sc_destroy_raster_env_includes(inc) ;
		sc_destroy_raster_bands_evaluator(eval) ;
	}
	lwpoly_free(outline) ;

	return sc ;
}

void
sc_destroy_raster_wrapper(SPATIAL_COLLECTION *dead)
{
	if (dead != NULL) {
		sc_destroy_raster_env_includes(dead->inclusion) ;
		sc_destroy_raster_bands_evaluator(dead->evaluator) ;
		sc_destroy(dead) ;
	}
}


SPATIAL_COLLECTION *
sc_create_raster_nodata_wrapper(rt_raster raster, int *bands, int num_bands,
		int nodata_band)
{
	INCLUDES *inc ;
	EVALUATOR *eval ;
	LWPOLY *outline ;
	SPATIAL_COLLECTION *sc ;

	if (raster==NULL) return NULL ;

	inc = sc_create_raster_nodata_includes(raster, nodata_band) ;
	if (inc == NULL) return NULL ;

	eval = sc_create_raster_bands_evaluator(raster, bands, num_bands) ;
	if (eval == NULL) {
		sc_destroy_raster_nodata_includes(inc) ;
		return NULL ;
	}

	outline = rt_raster_get_convex_hull(raster) ;
	lwgeom_add_bbox(lwpoly_as_lwgeom(outline)) ;

	sc = sc_create(SPATIAL_PLUS_VALUE,
			rt_raster_get_srid(raster), outline->bbox,
			NULL, inc, eval) ;

	if (sc == NULL) {
		sc_destroy_raster_nodata_includes(inc) ;
		sc_destroy_raster_bands_evaluator(eval) ;
	}
	lwpoly_free(outline) ;

	return sc ;
}

void
sc_destroy_raster_nodata_wrapper(SPATIAL_COLLECTION *dead)
{
	if (dead != NULL) {
		sc_destroy_raster_nodata_includes(dead->inclusion) ;
		sc_destroy_raster_bands_evaluator(dead->evaluator) ;
		sc_destroy(dead) ;
	}
}

/** @} */ /* end of documentation group spatial_collection_i */

