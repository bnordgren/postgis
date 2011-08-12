#ifndef RT_COLLECTION_H
#define RT_COLLECTION_H

/**
 * A method signature to be used as a "value callback" in conjunction with
 * rt_raster_raster_op_engine. A two_raster_value_op is responsible for
 * producing a value to store in the result raster given the value of the
 * cell in either or both of the input rasters. This is called for every
 * cell in the result raster (even those not included in the result).
 */
typedef void ((*two_raster_value_op)(double *,double*,double *,int,int,double)) ;

/**
 * A two_raster_value_op which populates the result raster with a mask
 * value (true/false). All cells
 * in the result raster have a value (none have a nodata value.)
 */
void
rt_raster_mask_value(double *r1_data, double *r2_data,
		double *result_data, int nbands, int included, double nodata);

/**
 * A two_raster_value_op which populates the result raster with pixel data
 * from r1 (if possible), or
 * r2 (if not). It is assumed that at least one of the input rasters
 * has data for the result pixel if it is included. If the cell is not
 * included in the result, the pixel is set to the nodata value. If there is
 * more than one band, all bands are assumed to share the same nodata value.
 */
void rt_raster_copy_first(double *r1_data, double *r2_data,
		double *result_data, int nbands, int included, double nodata) ;


/**
 * A method signature to be used in the calculation of the "spatial part"
 * of a two-raster spatial operation. This method accepts two boolean values,
 * each of which indicates whether the corresponding input raster contains
 * data over the output raster's cell. The return value indicates whether the
 * result cell is included in the result.
 */
typedef int ((*two_raster_spatial_op)(int,int)) ;

/**
 * A two_raster_spatial_op implementing the "intersection" spatial
 * operation.
 */
int rt_raster_raster_intersection(int r1, int r2)  ;

/**
 * A two_raster_spatial_op implementing the "union" spatial operation.
 */
int rt_raster_raster_union(int r1, int r2) ;

/**
 * A two_raster_spatial_op implementing the "difference" spatial operation.
 */
int rt_raster_raster_difference(int r1, int r2) ;

/**
 * A two_raster_spatial_op implementing the "symmetric difference"
 * spatial operation.
 */
int rt_raster_raster_symdifference(int r1, int r2) ;


/**
 * The two-raster engine which iterates over all the cells in the
 * result raster. Both the spatial operation function and the value
 * operation function should be provided.
 */
void rt_raster_raster_op_engine(
		rt_raster r1, char *srtext1,
		rt_raster r2, char *srtext2,
		rt_raster result, char *srtext_result,
		int mask,
		two_raster_spatial_op spatial_op,
		two_raster_value_op value_op) ;


/**
 * A method signature used for functions which calculate the extent
 * of a result raster. The inputs are the convex hulls of the two
 * input rasters, and the output is the size of the result. Both
 * input rasters should be in the same projection.
 */
typedef LWGEOM *((*rt_raster_grid_prep_op)(LWPOLY*,LWPOLY*)) ;


/**
 * Calculates the maximum size of the result grid using an
 * "intersection" operation.
 */
LWGEOM *rt_raster_prep_intersection_grid(LWPOLY *r1, LWPOLY *r2) ;

/**
 * Calculates the maximum size of the result grid using an
 * "difference" operation.
 */
LWGEOM *rt_raster_prep_difference_grid(LWPOLY *r1, LWPOLY *r2) ;

/**
 * Calculates the maximum size of the result grid using an
 * "symmetric difference" operation.
 */
LWGEOM *rt_raster_prep_symdifference_grid(LWPOLY *r1, LWPOLY *r2) ;

/**
 * Calculates the maximum size of the result grid using an
 * "union" operation.
 */
LWGEOM *rt_raster_prep_union_grid(LWPOLY *r1, LWPOLY *r2);

/**
 * Calculates the size of the resultant grid using the specified
 * operation. The extents of the two input rasters are projected into
 * the "destination" coordinate system prior to computing the result
 * grid, which is returned. The returned raster has geometric information
 * only: no bands and no pixel data.
 */
rt_raster rt_raster_raster_prep(
		rt_raster r1, projPJ r1_srs,
		rt_raster r2, projPJ r2_srs,
		int dest_srid, projPJ dest_srs,
		rt_raster_grid_prep_op grid_prep) ;
#endif
