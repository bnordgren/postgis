#include "liblwgeom.h"
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

/**
 * \defgroup rt_align_inc Wraps an Include interface, allowing access by raster indices.
 *
 * @{
 * Given a collection and a raster grid definition, this method wraps the
 * collection with includesIndex() and evaluateIndex() methods. Access to
 * the collection using the provided raster grid definition is then possible.
 * The provided raster may be empty, however, the raster and the collection
 * must have the same srid.
 */


struct rasterwrap_inc_s {
	INCLUDES *wrapped_inc ;
	rt_raster align_raster ;
};

/**
 * An implementation of the includes method which merely calls the
 * includes method of the wrapped object. This is reasonable because
 * we require that the raster's geocoordinates and the wrapped
 * includes have the same srid. (The collection makes sure of this.)
 */
int
collection_rasterwrap_includes(INCLUDES *inc, LWPOINT *point)
{
	struct rasterwrap_inc_s *params ;
	INCLUDES *wrapped ;

	if (inc == NULL) return 0 ;
	if (inc->params == NULL) return 0 ;

	params = (struct rasterwrap_inc_s *)(inc->params) ;
	wrapped = params->wrapped_inc ;
	if (wrapped == NULL || wrapped->includes == NULL) return 0 ;

	return wrapped->includes(wrapped, point) ;
}


/**
 * An implementation of includesIndex() which allows the user to
 * access the collection using the pixel indices of the supplied
 * raster.
 */
int
collection_rasterwrap_includes_idx(INCLUDES *inc, LWPOINT *point)
{
	struct rasterwrap_inc_s *params ;
	POINT2D   idx;
	POINT2D   geo ;
	LWPOINT   *geo_point ;
	int result ;

	if (inc == NULL || point==NULL) return 0 ;
	if (inc->params == NULL) return 0 ;
	if (inc->collection == NULL) return 0 ;

	params = (struct rasterwrap_inc_s *)(inc->params) ;
	if (params->align_raster == NULL) return 0;
	if (params->wrapped_inc == NULL) return 0 ;
	if (params->wrapped_inc->includes == NULL) return 0 ;

	/* compute geocoordinates using the raster's affine transform */
	lwpoint_getPoint2d_p(point, &idx) ;
	rt_raster_cell_to_geopoint(params->align_raster,
			idx.x, idx.y, &(geo.x), &(geo.y)) ;

	/* construct an LWPOINT with the geocoordinates */
	geo_point = lwpoint_make2d(inc->collection->srid, geo.x, geo.y) ;
	if (geo_point == NULL) return 0 ;

	/* call the includes method of the wrapped Includes interface */
	result = params->wrapped_inc->includes(params->wrapped_inc, geo_point) ;
	lwpoint_free(geo_point) ;

	return result ;
}

/**
 * Augments any Includes object with the ability to be accessed by
 * the indices of the specified raster. Note that the resultant object
 * returns whether the collection includes the specified raster indices,
 * not whether the raster includes the indices...
 *
 * @param inc the Includes object to augment
 * @param alignTo the raster containing the desired grid definition
 */
INCLUDES *
sc_create_rasterwrap_includes(INCLUDES *inc, rt_raster alignTo)
{
	struct rasterwrap_inc_s *params ;
	INCLUDES *wrapper ;

	if (inc == NULL || alignTo == NULL) return NULL ;

	params = (struct rasterwrap_inc_s *)rtalloc(sizeof(struct rasterwrap_inc_s)) ;
	if (params == NULL) return NULL ;

	params->wrapped_inc = inc ;
	params->align_raster = alignTo ;

	wrapper = inc_create(params,
			collection_rasterwrap_includes,
			collection_rasterwrap_includes_idx) ;

	if (wrapper == NULL) {
		rtdealloc(params) ;
		return NULL ;
	}

	return wrapper ;
}

void sc_destroy_rasterwrap_includes(INCLUDES *dead)
{
	if (dead != NULL) {
		if (dead->params != NULL) {
			rtdealloc(dead->params) ;
		}
		inc_destroy(dead) ;
	}
}

/** @} */ /* end of documentation group rt_align_inc */
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

/**
 * \defgroup rt_align_eval Provides access to a collection using a raster's indicies.
 *
 * @{
 */

struct rasterwrap_eval_s {
	EVALUATOR *wrapped_eval ;
	rt_raster align_raster ;
};

/**
 * An implementation of the evaluate method of the Evaluator interface
 * which simply calls the evaluate method of the wrapped object. This requires
 * that the srid of this object, the srid of the wrapped object, and the srid
 * of the raster's geocoordinates all be identical.
 */
VALUE *
collection_rasterwrap_evaluate(EVALUATOR *eval, LWPOINT *point)
{
	struct rasterwrap_eval_s *params ;
	EVALUATOR *wrapped ;
	VALUE *wrapped_result ;

	if (eval == NULL ) return NULL ;
	if (eval->params == NULL) return NULL ;

	params = (struct rasterwrap_eval_s *)(eval->params) ;
	wrapped = params->wrapped_eval ;

	if (wrapped == NULL || wrapped->evaluate == NULL) return NULL ;

	wrapped_result = wrapped->evaluate(wrapped, point) ;
	val_copy(eval->result, wrapped_result) ;

	return eval->result ;
}

/**
 * An implementation of evaluateIndex() which allows the user to
 * access the collection using the pixel indices of the supplied
 * raster.
 */
VALUE *
collection_rasterwrap_evaluate_idx(EVALUATOR *eval, LWPOINT *point)
{
	struct rasterwrap_eval_s *params ;
	POINT2D   idx;
	POINT2D   geo ;
	LWPOINT   *geo_point ;
	VALUE     *wrapped_result ;

	if (eval == NULL || point==NULL) return NULL ;
	if (eval->params == NULL) return NULL ;
	if (eval->collection == NULL) return NULL ;

	params = (struct rasterwrap_eval_s *)(eval->params) ;
	if (params->align_raster == NULL) return NULL;
	if (params->wrapped_eval == NULL) return NULL ;
	if (params->wrapped_eval->evaluate == NULL) return NULL ;

	/* compute geocoordinates using the raster's affine transform */
	lwpoint_getPoint2d_p(point, &idx) ;
	rt_raster_cell_to_geopoint(params->align_raster,
			idx.x, idx.y, &(geo.x), &(geo.y)) ;

	/* construct an LWPOINT with the geocoordinates */
	geo_point = lwpoint_make2d(eval->collection->srid, geo.x, geo.y) ;
	if (geo_point == NULL) return NULL ;

	/* call the evaluate method of the wrapped Evaluator interface */
	wrapped_result = params->wrapped_eval->evaluate(params->wrapped_eval, geo_point) ;
	val_copy(eval->result, wrapped_result) ;
	lwpoint_free(geo_point) ;

	return eval->result ;
}

EVALUATOR *
sc_create_rasterwrap_evaluator(EVALUATOR *eval, rt_raster alignTo)
{
	struct rasterwrap_eval_s *params ;
	EVALUATOR *wrapper ;

	if (eval == NULL || alignTo == NULL) return NULL ;
	if (eval->result == NULL) return NULL ;


	params = (struct rasterwrap_eval_s *)rtalloc(sizeof(struct rasterwrap_eval_s)) ;
	if (params == NULL) return NULL ;

	params->wrapped_eval = eval ;
	params->align_raster = alignTo ;

	wrapper = eval_create(params,
			collection_rasterwrap_evaluate,
			collection_rasterwrap_evaluate_idx,
			eval->result->length) ;

	if (wrapper == NULL) {
		rtdealloc(params) ;
		return NULL ;
	}

	return wrapper ;
}

void
sc_destroy_rasterwrap_evaluator(EVALUATOR *dead)
{
	if (dead != NULL) {
		if (dead->params != NULL) {
			rtdealloc(dead->params) ;
		}
		eval_destroy(dead) ;
	}
}

/** @} */ /* end of documentation group rt_align_eval */

/** @} */ /* end of documentation group evaluator_i */

/**
 * \addtogroup spatial_collection_i Implementations of the "SpatialCollection" interface.
 * @{
 */


/**
 * \defgroup sc_rasterwrap Raster Wrapper
 *
 * @{
 * An implementation of #SPATIAL_COLLECTION which is backed by a raster.
 * The raster is assumed to lack a "nodata" value, which means that its
 * spatial extent is identical to its convex hull. (e.g., every grid cell
 * has a value). The value returned by the collection is the value of one or
 * more of the bands at the evaluation point.
 */

struct raster_wrap_s {
	rt_raster raster ;
	int       owned ;
};

/**
 * Constructs a #SPATIAL_COLLECTION backed by a raster without a nodata
 * value. The evaluate method on this collection returns the specified
 * bands.
 *
 * The collection returned by this constructor should be freed by
 * #sc_destroy_raster_wrapper.
 *
 * @param raster the raster to wrap
 * @param owned  if TRUE, the wrapper takes ownership of the raster and
 *               will free it when the wrapper is destroyed
 * @param bands  a list of the bands, in order, which should be returned
 *               from the evaluate method.
 * @param num_bands the length of bands
 */
SPATIAL_COLLECTION *
sc_create_raster_wrapper(rt_raster raster, int owned, int *bands, int num_bands)
{
	INCLUDES *inc ;
	EVALUATOR *eval ;
	LWPOLY *outline ;
	SPATIAL_COLLECTION *sc ;
	struct raster_wrap_s *params ;

	if (raster==NULL) return NULL ;

	inc = sc_create_raster_env_includes(raster) ;
	if (inc == NULL) return NULL ;

	eval = sc_create_raster_bands_evaluator(raster, bands, num_bands) ;
	if (eval == NULL) {
		sc_destroy_raster_env_includes(inc) ;
		return NULL ;
	}

	params = (struct raster_wrap_s *)rtalloc(sizeof(struct raster_wrap_s));
	if (params == NULL) {
		sc_destroy_raster_env_includes(inc) ;
		sc_destroy_raster_bands_evaluator(eval) ;
		return NULL ;
	}
	params->raster = raster ;
	params->owned  = owned ;

	outline = rt_raster_get_convex_hull(raster) ;
	lwgeom_add_bbox(lwpoly_as_lwgeom(outline)) ;

	sc = sc_create(SPATIAL_PLUS_VALUE,
			rt_raster_get_srid(raster),
			outline->bbox,
			params, inc, eval) ;

	if (sc == NULL) {
		sc_destroy_raster_env_includes(inc) ;
		sc_destroy_raster_bands_evaluator(eval) ;
	}
	lwpoly_free(outline) ;

	return sc ;
}

/**
 * Destroys a #SPATIAL_COLLECTION object which was constructed by
 * #sc_create_raster_wrapper or #sc_create_pgraster_wrapper.
 * Pass no other objects to this function.
 */
void
sc_destroy_raster_wrapper(SPATIAL_COLLECTION *dead)
{
	struct raster_wrap_s *params ;
	if (dead != NULL) {
		sc_destroy_raster_env_includes(dead->inclusion) ;
		sc_destroy_raster_bands_evaluator(dead->evaluator) ;

		if (dead->params != NULL) {
			params = (struct raster_wrap_s *)(dead->params) ;
			if ((params->raster != NULL) && params->owned) {
				rt_raster_destroy(params->raster) ;
			}
			rtdealloc(params) ;
		}

		sc_destroy(dead) ;
	}
}

/** @} */ /* end of sc_rasterwrap documentation group */

/**
 * \defgroup sc_rasterwrap_nodata Raster Wrapper w/Nodata
 *
 * @{
 */

/**
 * Constructs a #SPATIAL_COLLECTION backed by a raster with a nodata
 * value. The evaluate method on this collection returns the specified
 * bands.
 *
 * The collection returned by this constructor should be freed by
 * #sc_destroy_raster_nodata_wrapper.
 *
 * @param raster the raster to wrap
 * @param owned  if TRUE, the wrapper takes ownership of the raster and
 *               will free it when the wrapper is destroyed
 * @param bands  a list of the bands, in order, which should be returned
 *               from the evaluate method.
 * @param num_bands the length of bands
 * @param nodata_band the index of the band which will be used to test
 *               for a nodata value.
 */
SPATIAL_COLLECTION *
sc_create_raster_nodata_wrapper(rt_raster raster, int owned,
		int *bands, int num_bands, int nodata_band)
{
	INCLUDES *inc ;
	EVALUATOR *eval ;
	LWPOLY *outline ;
	SPATIAL_COLLECTION *sc ;
	struct raster_wrap_s *params ;

	if (raster==NULL) return NULL ;

	inc = sc_create_raster_nodata_includes(raster, nodata_band) ;
	if (inc == NULL) return NULL ;

	eval = sc_create_raster_bands_evaluator(raster, bands, num_bands) ;
	if (eval == NULL) {
		sc_destroy_raster_nodata_includes(inc) ;
		return NULL ;
	}

	params = (struct raster_wrap_s *)rtalloc(sizeof(struct raster_wrap_s)) ;
	if (params == NULL) {
		sc_destroy_raster_nodata_includes(inc) ;
		sc_destroy_raster_bands_evaluator(eval) ;
		return NULL ;
	}
	params->raster = raster ;
	params->owned  = owned ;

	outline = rt_raster_get_convex_hull(raster) ;
	lwgeom_add_bbox(lwpoly_as_lwgeom(outline)) ;

	sc = sc_create(SPATIAL_PLUS_VALUE,
			rt_raster_get_srid(raster), outline->bbox,
			params, inc, eval) ;

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

/** @} */ /* end of documentation group sc_rasterwrap_nodata */

/**
 * \defgroup collection_rasterwrap Spatial collection providing access via the indices of the provided raster.
 *
 * @{
 * This spatial collection passes the includes() and evaluates() calls thru to the
 * wrapped collection. The includesIndex() and evaluateIndex() calls use coordinates
 * aligned to the grid defined by the provided raster.
 */

SPATIAL_COLLECTION *
sc_create_raster_aligned_collection(SPATIAL_COLLECTION *wrapped, rt_raster alignTo)
{
	EVALUATOR *eval ;
	INCLUDES *inc ;
	SPATIAL_COLLECTION *wrapper ;
	if (wrapped == NULL || alignTo==NULL) return NULL  ;
	if (wrapped->inclusion == NULL) return NULL ;

	inc = sc_create_rasterwrap_includes(wrapped->inclusion, alignTo) ;
	if (inc == NULL) return NULL ;

	eval = NULL  ;
	if (wrapped->evaluator != NULL) {
		eval = sc_create_rasterwrap_evaluator(wrapped->evaluator, alignTo) ;
	}

	wrapper = sc_create(wrapped->type,
			wrapped->srid,
			&(wrapped->extent),
			NULL, inc, eval) ;

	if (wrapper == NULL) {
		sc_destroy_rasterwrap_includes(inc) ;
		if (eval != NULL) {
			sc_destroy_rasterwrap_evaluator(eval) ;
		}
		return NULL ;
	}

	return wrapper ;

}

void
sc_destroy_raster_aligned_collection(SPATIAL_COLLECTION *dead)
{
	if (dead != NULL) {
		sc_destroy_rasterwrap_includes(dead->inclusion) ;
		if (dead->evaluator != NULL) {
			sc_destroy_rasterwrap_evaluator(dead->evaluator) ;
		}
		sc_destroy(dead) ;
	}
}

/** @} */  /* end of collection_rasterwrap documentation group */


/** @} */ /* end of documentation group spatial_collection_i */

/**
 * Samples a #SPATIAL_COLLECTION at each grid cell in the provided
 * raster. The raster must have a grid definition but if it lacks bands,
 * the bands will be created. If the raster does have bands, the number of
 * bands must be the same as the length of the #VALUE object returned by the
 * evaluate method on the spatial collection.
 *
 * The user may provide a nodata vector to be written in each cell of
 * the raster not included in the result. If the user does not provide
 * a nodata #VALUE, then zero will be written to all bands.
 *
 * @param source The spatial collection to be sampled
 * @param result The raster containing grid points at which the sampling
 *               takes place.
 * @param nodata_val (optional) The value to write to a grid cell if
 *                   the cell is not included in the result.
 */
void
sc_sampling_engine(SPATIAL_COLLECTION *source,
		           rt_raster result,
		           VALUE *nodata_val)
{
	uint16_t width ;
	uint16_t height ;
	VALUE   *collection_val ;
	int      coll_bands ;
	int      raster_bands ;
	int      local_nodata ;
	int      i,j ;
	LWPOINT *sample_pt ;
	POINT4D  sample_pt_p4d ;
	SPATIAL_COLLECTION *aligned ;

	if (source == NULL || result == NULL) return ;
	if (source->type == SPATIAL_ONLY) {
		rterror("sc_sampling_engine: cannot sample a collection with no values.");
		return ;
	}
	if (source->evaluator == NULL) return ;
	if (source->evaluator->result == NULL) return ;

	/* get dimensions of the result. */
	width = result->width ;
	height = result->height;

	/* align the collection to the result raster */
	aligned = sc_create_raster_aligned_collection(source, result);
	if (aligned == NULL) return ;

	/* get number of bands */
	collection_val = aligned->evaluator->result ;
	coll_bands = collection_val->length ;
	raster_bands = rt_raster_get_num_bands(result) ;

	/* is the raster ready to receive values from the collection? */
	if ((raster_bands != 0) && (raster_bands != coll_bands)) {
		rterror("sc_sampling_engine: non-empty raster must have same number of bands as the sampled collection.") ;
		return ;
	}

	/* did the user specify a nodata value? */
	local_nodata = (nodata_val == NULL) ;
	if (local_nodata) {
		int nd_band ;

		nodata_val = val_create(coll_bands) ;
		rterror("sc_sampling_engine: cannot create NODATA vector and none specified") ;
		if (nodata_val == NULL) return ;

		/* initialize nodata vector to all zeros */
		for (nd_band=0; nd_band < coll_bands; nd_band++) {
			nodata_val->data[nd_band] = 0.0 ;
		}
	}

	/* maybe we have to add the bands ourselves */
	if (raster_bands == 0) {
		int band ;
		for (band=0; band < coll_bands; band++) {
			uint8_t *band_data ;
			rt_band empty ;
			int     pix_size ;
			rt_pixtype band_type ;

			/* note: need to make VALUE keep track of datatype */
			band_type = PT_64BF ;
			pix_size = rt_pixtype_size(band_type) ;
			band_data = rtalloc(pix_size * width * height) ;
			if (band_data == NULL) {
				rterror("sc_sampling_engine: cannot allocate memory for band") ;

				/* free up successfully allocated bands */
				band -- ;
				while (band >= 0) {
					empty = rt_raster_get_band(result, band) ;
					band_data = (uint8_t *)rt_band_get_data(empty);
					rt_band_destroy(empty) ;
					band -- ;
				}
				return ;
			}

			/* make the new band and add it to the raster. */
			empty = rt_band_new_inline(width, height, band_type, 0, 0.0, band_data);
			rt_raster_add_band(result, empty, band) ;
		}
	}

	/* now sample the raster */
	sample_pt = lwpoint_make2d(SRID_UNKNOWN, 0,0) ;
	for (j=0; j<height; j++) {
		for (i=0; i<width; i++) {
			int band_num ;
			VALUE *store_val ;

			/* set up the sample point */
			sample_pt_p4d.x = i ;
			sample_pt_p4d.y = j ;
			ptarray_set_point4d(sample_pt->point,0,&sample_pt_p4d) ;

			/* evaluate the collection */
			collection_val = sc_evaluateIndex(aligned, sample_pt) ;

			/* store the returned value or the nodata vector */
			store_val = collection_val ;
			if (collection_val == NULL) {
				store_val = nodata_val ;
			}

			/* copy the values to the raster */
			for (band_num=0; band_num < coll_bands; band_num++) {
				rt_band band ;

				band = rt_raster_get_band(result, band_num) ;
				rt_band_set_pixel(band, i, j, store_val->data[band_num]) ;
			}
		}
	}

	sc_destroy_raster_aligned_collection(aligned) ;
	if (local_nodata) {
		val_destroy(nodata_val) ;
	}

}

/**
 * Assuming the presence of a computed result extent (as in the envelope
 * of the result of a spatial operation) and a raster having the four
 * orientation parameters set, this function calculates the offsets required
 * to place the raster in the desired extent. The method leverages the
 * fact that the size of the raster's extent will be the same regardless of
 * where it is. It calculates the xmin, ymin corner of the extent when
 * the offsets are set to zero (at the origin). Then it translates the
 * whole raster (including the origin) by the amount necessary to put the
 * calculated xmin, ymin corner on top of the xmin, ymin corner of the
 * provided extent.
 *
 * This method works for any rotation angle, but requires that the
 * width, height, and scales of the raster be consistent with the
 * dimensions of the extent.
 *
 * @param extent the extent in which the raster should live
 * @param raster the raster which needs to be translated into the extent
 */
void
fit_raster_to_extent(GBOX *extent, rt_raster raster)
{
	double o11, o12, o21, o22 ;
	double x[3] ;
	double y[3] ;
	double xmin, ymin ;
	double xoffset, yoffset ;
	int width, height ;
	int i;

	if (extent==NULL || raster==NULL) return ;

	width = rt_raster_get_width(raster);
	height = rt_raster_get_height(raster) ;

	/* get the matrix coefficients */
	o11 = rt_raster_get_scaleX(raster) ;
	o12 = rt_raster_get_skewX(raster) ;
	o21 = rt_raster_get_skewY(raster) ;
	o22 = rt_raster_get_scaleY(raster) ;

	/* upper-left corner is origin (0,0) */

	/* upper-right corner (pixel (width,0)) */
	x[0] = o11 * width ;
	y[0] = o21 * width ;

	/* lower-right corner (pixel (width, height)) */
	x[1] = o11 * width + o12 * height ;
	y[1] = o21 * width + o22 * height ;

	/* lower left corner (pixel (0, height)) */
	x[2] = o12 * height ;
	y[2] = o22 * height ;

	xmin = ymin = 0 ; /* origin */
	for (i=0; i<3 ; i++) {
		xmin = fmin(xmin, x[i]) ;
		ymin = fmin(ymin, y[i]) ;
	}

	/* compute the offset of the "extent" box centered around 0,0 */
	xoffset = extent->xmin - xmin ;
	yoffset = extent->ymin - ymin ;

	rt_raster_set_offsets(raster, xoffset, yoffset) ;
}
