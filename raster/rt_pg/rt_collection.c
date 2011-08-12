#include "liblwgeom.h"
#include "rt_api.h"
#include "rt_core_internal.h"
#include "rt_collection.h"

/*- spatial operations (spatial part) --------------------------------*/

int rt_raster_raster_intersection(int r1, int r2)
{
	return r1 && r2 ;
}

int rt_raster_raster_union(int r1, int r2)
{
	return r1 || r2 ;
}

int rt_raster_raster_difference(int r1, int r2)
{
	return r1 && !r2  ;
}

int rt_raster_raster_symdifference(int r1, int r2)
{
	return (r1 && !r2) || (!r1 && r2) ;
}

/*- spatial operations (value part) ----------------------------------*/

void rt_raster_mask_value(double *r1_data, double *r2_data,
		double *result_data, int nbands, int included, double nodata) {
	result_data[0] = included ;
}

void rt_raster_copy_first(double *r1_data, double *r2_data,
		double *result_data, int nbands, int included, double nodata) {
	int band ;
	for (band = 0; band < nbands ; band ++) {
		if (!included || (!r1_data && !r2_data)) {
			result_data[band] = nodata ;
		} else if (r1_data) {
			result_data[band] = r1_data[band] ;
		} else if (r2_data) {
			result_data[band] = r2_data[band] ;
		}
	}
}

/*- spatial operations (raster iteration) -----------------------------*/
void rt_raster_raster_op_engine(
		rt_raster r1, char *srtext1,
		rt_raster r2, char *srtext2,
		rt_raster result, char *srtext_result,
		int mask,
		two_raster_spatial_op spatial_op,
		two_raster_value_op value_op) {

	int band ; /* for iteration */

	/* get dimensions of the result. */
	uint16_t width = result->width ;
	uint16_t height = result->height;

	/* get band counts for input rasters */
	int num_bands_r1, num_bands_r2, num_bands_result;
	if (mask) {
		/* to generate mask, only need first band */
		num_bands_r1 = 1 ;
		num_bands_r2 = 1 ;
		num_bands_result = 1 ;
	} else {
		num_bands_r1 = rt_raster_get_num_bands(r1) ;
		num_bands_r2 = rt_raster_get_num_bands(r2) ;
		num_bands_result = rt_raster_get_num_bands(result) ;
	}

	/* check that the number of bands is the same for all rasters */
	if ( (num_bands_r1 != num_bands_r2) ||
			(num_bands_r2 != num_bands_result)) {
		rterror("rt_raster_raster_op_engine: all rasters must have same number of bands") ;
		return ;
	}
	/* check that rasters have at least one band */
	if (num_bands_r1 <=0) {
		rterror("rt_raster_raster_op_engine: rasters have no band data");
		return ;
	}


	/* get the bands for the rasters */
	rt_band *r1_bands, *r2_bands, *result_bands ;
	r1_bands = (rt_band *)(rtalloc(sizeof(rt_band)*num_bands_r1*3)) ;
	r2_bands = r1_bands + num_bands_r1 ;
	result_bands = r2_bands + num_bands_r2 ;
	for (band=0; band<num_bands_r1; band++) {
		r1_bands[band] = rt_raster_get_band(r1, band) ;
		r2_bands[band] = rt_raster_get_band(r2, band) ;
		result_bands[band] = rt_raster_get_band(result, band) ;
	}

	/* check that result raster has correct nodata */
	if (!rt_band_get_hasnodata_flag(result_bands[0]) && !mask) {
		rterror("rt_raster_raster_op_engine: non-mask result must have nodata value") ;
		rtdealloc(r1_bands) ;
		return ;
	}
	if (rt_band_get_hasnodata_flag(result_bands[0]) && mask) {
		rterror("rt_raster_raster_op_engine: mask result must not have a nodata value");
		rtdealloc(r1_bands) ;
		return ;
	}
	double result_nodata_val = 0.;
	if (!mask) result_nodata_val = rt_band_get_nodata(result_bands[0]) ;

	/* get nodata information for r1 and r2 */
	int r1_has_nodata = rt_band_get_hasnodata_flag(r1_bands[0]) ;
	int r2_has_nodata = rt_band_get_hasnodata_flag(r2_bands[0]) ;
	double r1_nodata_val,r2_nodata_val = 0. ;
	if (r1_has_nodata) {
		r1_nodata_val = rt_band_get_nodata(r1_bands[0]) ;
	}
	if (r2_has_nodata) {
		r2_nodata_val = rt_band_get_nodata(r2_bands[0]) ;
	}


	/* create buffers for single pixel (all bands) */
	double *r1_data, *r2_data, *result_data ;
	r1_data = (double *)rtalloc(sizeof(double)*num_bands_r1*3) ;
	r2_data = r1_data + num_bands_r1 ;
	result_data = r2_data + num_bands_r1 ;

	/* Get GDALDatasets for inputs and result */
	GDALDriverH drv ;
	uint32_t firstband = 0 ;
	GDALDatasetH gdal_r1 = rt_raster_to_gdal_mem(r1, srtext1,
			&firstband, 1, &drv) ;
	GDALDatasetH gdal_r2 = rt_raster_to_gdal_mem(r2, srtext2,
			&firstband, 1, &drv) ;
	GDALDatasetH gdal_result = rt_raster_to_gdal_mem(result, srtext_result,
			&firstband, 1, &drv) ;

	/* Create transformation between result and the two inputs. */
	void *r1_xform = GDALCreateGenImgProjTransformer2(
			gdal_result, gdal_r1, NULL) ;
	void *r2_xform = GDALCreateGenImgProjTransformer2(
			gdal_result, gdal_r2, NULL) ;

	double i_r1, j_r1, i_r2, j_r2, k ;
	int    r1_good, r2_good ;
	int    i,j ;
	for (i=0; i < width; i ++) {
		for (j=0; j<height; j++) {
			i_r1 = i_r2 = i ;
			j_r1 = j_r2 = j ;
			k = 0.0 ;
			r1_good = r2_good = 0 ;

			/* calculate r1 image coordinates */
			GDALGenImgProjTransform(r1_xform, 0 /* false: src->dst */,
					1, &i_r1, &j_r1, &k, &r1_good) ;

			/* calculate r2 image coordinates */
			k = 0.0 ;
			GDALGenImgProjTransform(r2_xform, 0 /* false: src->dst */,
					1, &i_r2, &j_r2, &k, &r2_good) ;

			/* get input data */
			uint16_t tmp_i_r1 = (uint16_t)(round(i_r1)) ;
			uint16_t tmp_j_r1 = (uint16_t)(round(j_r1)) ;
			uint16_t tmp_i_r2 = (uint16_t)(round(i_r2)) ;
			uint16_t tmp_j_r2 = (uint16_t)(round(j_r2)) ;
			int r1_inrange = r1_good &&
					((tmp_i_r1 >= 0) && (tmp_i_r1 < r1->width)) &&
					((tmp_j_r1 >= 0) && (tmp_j_r1 < r1->height)) ;
			int r2_inrange = r2_good &&
					((tmp_i_r2 >= 0) && (tmp_i_r2 < r2->width)) &&
					((tmp_j_r2 >= 0) && (tmp_j_r2 < r1->height)) ;

			for (band = 0; band < num_bands_r1; band++) {
				if (r1_inrange) {
					rt_band_get_pixel(r1, tmp_i_r1, tmp_j_r1, r1_data+band);
				}
				if (r2_inrange) {
					rt_band_get_pixel(r2, tmp_i_r2, tmp_j_r2, r2_data+band);
				}
			}

			/* compute whether the input rasters contain data
			 * over the "result" cell.
			 */
			int r1_present = r1_inrange &&
				(!r1_has_nodata ||
				 (r1_has_nodata && (r1_data[0]!=r1_nodata_val)));
			int r2_present = r2_inrange &&
				(!r2_has_nodata ||
				 (r2_has_nodata && (r2_data[0]!=r2_nodata_val)));

			/* compute whether the result cell is inside or outside the
			 * spatial result.
			 */
			int included = spatial_op(r1_present, r2_present) ;

			/* now act on the result of the spatial op */
			value_op( (r1_present) ? r1_data : NULL,
					  (r2_present) ? r2_data : NULL,
					  result_data,
					  num_bands_result, included, result_nodata_val) ;

			/* write the result to the output */
			for (band=0; band < num_bands_result; band++) {
				rt_band_set_pixel(result_bands[band], i, j,
						result_data[band]);
			}
		}
	}


	/* Free memory */
	GDALDestroyGenImgProjTransformer(r1_xform) ;
	GDALDestroyGenImgProjTransformer(r2_xform) ;
	GDALClose(gdal_r1) ;
	GDALClose(gdal_r2) ;
	GDALClose(gdal_result) ;
	rtdealloc(r1_bands) ;
	rtdealloc(r1_data) ;

}

/*- spatial operations (raster preparation) --------------------------*/

LWGEOM *rt_raster_prep_intersection_grid(LWPOLY *r1, LWPOLY *r2)
{
	return lwgeom_intersection(lwpoly_as_lwgeom(r1),
			                   lwpoly_as_lwgeom(r2)) ;
}

LWGEOM *rt_raster_prep_difference_grid(LWPOLY *r1, LWPOLY *r2)
{
	return lwgeom_difference(lwpoly_as_lwgeom(r1),
			                 lwpoly_as_lwgeom(r2)) ;
}

LWGEOM *rt_raster_prep_union_grid(LWPOLY *r1, LWPOLY *r2)
{
	return lwgeom_union(lwpoly_as_lwgeom(r1),
			            lwpoly_as_lwgeom(r2)) ;
}

LWGEOM *rt_raster_prep_symdifference_grid(LWPOLY *r1, LWPOLY *r2)
{
	return lwgeom_symdifference(lwpoly_as_lwgeom(r1),
			                    lwpoly_as_lwgeom(r2)) ;
}

rt_raster rt_raster_raster_prep(
		rt_raster r1, projPJ r1_srs,
		rt_raster r2, projPJ r2_srs,
		int dest_srid, projPJ dest_srs,
		rt_raster_grid_prep_op grid_prep)
{

	/* get basic spatial info about raster extents */
	int r1_srid = rt_raster_get_srid(r1) ;
	int r2_srid = rt_raster_get_srid(r2) ;
	LWPOLY *r1_outline = rt_raster_get_convex_hull(r1);
	LWPOLY *r2_outline = rt_raster_get_convex_hull(r2);


	/* project outline of r1 if necessary */
	if (r1_srid != dest_srid) {
		lwgeom_transform(r1_outline, r1_srs, dest_srs) ;
	}

	/* project outline of r2 if necessary */
	if (r2_srid != dest_srid) {
		lwgeom_transform(r2_outline, r2_srs, dest_srs) ;
	}

	/* compute the extent of the result */
	LWGEOM *result_outline = grid_prep(r1_outline, r2_outline) ;
	lwgeom_add_bbox(result_outline) ;

	/* construct the raster to hold the result */
	GBOX *result_extent = result_outline->bbox ;
	double xres = fmin(rt_raster_get_x_scale(r1), rt_raster_get_x_scale(r2)) ;
	double yres = fmin(rt_raster_get_y_scale(r1), rt_raster_get_y_scale(r2)) ;
	uint16_t width = (int)ceil( (result_extent->xmax-result_extent->xmin)/xres) ;
	uint16_t height = (int)ceil( (result_extent->ymax-result_extent->ymin)/yres) ;
	rt_raster result = rt_raster_new(width,height) ;
	rt_raster_set_offsets(result, result_extent->xmin, result_extent->ymax) ;
	rt_raster_set_scale(result, xres, yres) ;
	result->srid = dest_srid ;

	lwgeom_free(r1_outline) ;
	lwgeom_free(r2_outline) ;
	lwgeom_free(result_outline) ;

	return result ;
}
