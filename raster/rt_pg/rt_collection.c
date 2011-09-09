#include "liblwgeom.h"
#include "rt_api.h"
#include "spatial_collection.h"
#include "lwgeom_pg.h"
#include "rt_collection.h"

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
	LWGEOM geom ;

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

