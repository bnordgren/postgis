#include "liblwgeom.h"
#include "rt_api.h"
#include "spatial_collection.h"
#include "lwgeom_pg.h"
#include "sc_raster.h"
#include "rt_collection.h"
#include "postgres.h"
#include "fmgr.h"
#include "utils/array.h"

/**
 * \addtogroup geo_wrap_collection pglwgeom wrapper
 * @{
 */

/**
 * This is an additional constructor of a geometry wrapper collection
 * object. Instead of accepting an #LWGEOM, it accepts a #PG_LWGEOM,
 * deserializes it, and wraps the resultant #LWGEOM object (which is
 * hidden from the caller.)
 *
 * The spatial collection returned from this constructor should be
 * destroyed with @sc_destroy_geometry_wrapper.
 */
SPATIAL_COLLECTION *
sc_create_pglwgeom_wrapper(PG_LWGEOM *pg_geom, double inside, double outside)
{
	LWGEOM *geom ;

	if (pg_geom == NULL) return NULL ;
	geom = pglwgeom_deserialize(pg_geom) ;
	return sc_create_geometry_wrapper(geom, 1, inside, outside) ;
}

/** @} */ /* end of geo_wrap_collection documentation group */

/**
 * \addtogroup sc_rasterwrap
 * @{
 */

/**
 * This is an additional constructor for the raster wrapper
 * spatial collection implementation. This constructor takes the
 * serialized form of the raster, deserializes it, and backs the
 * collection with the deserialized version.
 *
 * The object returned by this constructor should be destroyed using
 * #sc_destroy_raster_wrapper.
 *
 * @params pg_raster the serialized raster
 * @params bands     list of bands containing values which will be returned
 *                   by the evaluate method.
 * @params num_bands length of bands
 *
 */
SPATIAL_COLLECTION *
sc_create_pgraster_wrapper(rt_pgraster *pg_raster, int *bands, int num_bands)
{
	rt_raster raster ;

	raster = rt_raster_deserialize(pg_raster, 0) ;
	return sc_create_raster_wrapper(raster, 1, bands, num_bands) ;
}
/** @} */ /* end of sc_rasterwrap documentation group */

/**
 * \addtogroup sc_rasterwrap_nodata
 * @{
 */

/**
 * This is an additional constructor for the raster wrapper
 * spatial collection implementation, where the raster is assumed to
 * possess a nodata value. This constructor takes the
 * serialized form of the raster, deserializes it, and backs the
 * collection with the deserialized version.
 *
 * The object returned by this constructor should be destroyed using
 * #sc_destroy_raster_nodata_wrapper.
 *
 * @params pg_raster the serialized raster
 * @params bands     list of bands containing values which will be returned
 *                   by the evaluate method.
 * @params num_bands length of bands
 * @param nodata_band the index of the band which will be used to test
 *               for a nodata value.
 *
 */
SPATIAL_COLLECTION *
sc_create_pgraster_wrapper_nodata(rt_pgraster *pg_raster,
		int *bands, int num_bands, int nodata_band)
{
	rt_raster raster ;

	raster = rt_raster_deserialize(pg_raster, 0) ;
	return sc_create_raster_nodata_wrapper(raster, 1,
			bands, num_bands, nodata_band) ;
}

/** @} */ /* end of sc_rasterwrap_nodata documentation group */


/**
 * This function takes an array from argnum position of the function call
 * list and sets the pointer to the array of bands as well as the
 * number of bands. Consider this read only. Free the memory
 * returned in bands using rtdealloc. This buffer is part of the FunctionCallInfo...
 *
 * If the array is NULL, this function returns a list of all the
 * bands in the provided raster.
 */
int
getarg_bandlist(FunctionCallInfo fcinfo, int argnum,
		        rt_pgraster *raster, int **bands, int *num_bands)
{

	if (bands == NULL || num_bands==NULL || raster==NULL) return 0 ;
	if (PG_ARGISNULL(argnum)) {
		int i ;

		*num_bands = raster->numBands ;
		*bands = (int*)rtalloc(sizeof(int) * (*num_bands)) ;
		if (*bands == NULL) {
			rterror("getarg_bandlist: cannot allocate default band list.") ;
			return 0 ;
		}
		for (i=0; i<(*num_bands) ; i++) {
			(*bands)[i] = i ;
		}

	} else {
		ArrayType *pg_bands ;

		/* get the array and check for null */
		pg_bands = PG_GETARG_ARRAYTYPE_P(argnum) ;
		if (ARR_HASNULL(pg_bands)) {
			rterror("getarg_bandlist: NULL values not allowed in band list.") ;
			return 0 ;
		}

		/* check that the array is a vector */
		if (ARR_NDIM(pg_bands) != 1) {
			rterror("getarg_bandlist: band list array must be one-dimensional.") ;
			return 0 ;
		}

		/* check that the array contains integers */
		// dunno how to do this
		// use ARR_ELEMTYPE(), but how to divine the Oid associated with
		// integers?

		/* get the number of bands in the band list */
		*num_bands = (ARR_DIMS(pg_bands))[0] ;

		/* allocate the proper amount of memory in the band buffer */
		*bands = (int *)rtalloc(sizeof(int) * (*num_bands)) ;
		memcpy(*bands, ARR_DATA_PTR(pg_bands), sizeof(int)* (*num_bands)) ;
	}

	return 1;
}

/**
 * Creates a new, empty raster with the same rotation and scale
 * present in the grid_defn.
 *
 * @param extent the area which should be covered by the new raster
 *               (e.g., a box which should be inscribed in the raster)
 * @param grid_defn a grid from which we take alignment information.
 * @returns a raster of the correct size, or null if an error.
 */
rt_raster
rt_raster_new_inbox(GBOX *extent, rt_pgraster *grid_defn)
{
	uint16_t res_width, res_height ;
	rt_raster result ;

	if (extent == NULL || grid_defn == NULL) return NULL ;

	/* I know this is wrong for the moment. Need to get real calculation */
	res_width = (int)fabs((extent->xmax-extent->xmin) / grid_defn->scaleX) ;
	res_height = (int)fabs((extent->ymax-extent->ymin) / grid_defn->scaleY) ;

	/* make an empty raster to store the result */
	result = rt_raster_new(res_width, res_height) ;

	/* copy srid and geo transform (except for offsets) */
	rt_raster_set_srid(result, grid_defn->srid) ;
	rt_raster_set_scale(result, grid_defn->scaleX, grid_defn->scaleY) ;
	rt_raster_set_skews(result, grid_defn->skewX, grid_defn->skewY) ;

	/* compute new offsets because raster may be rotated
	 * w.r.t. the extent.
	 */
	fit_raster_to_extent(extent, result) ;

	return result ;
}


