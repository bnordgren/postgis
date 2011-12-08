#ifndef SC_RASTER_H
#define SC_RASTER_H

#include "liblwgeom.h"
#include "spatial_collection.h"
#include "rt_api.h"

INCLUDES *
sc_create_raster_env_includes(rt_raster raster) ;
void
sc_destroy_raster_env_includes(INCLUDES *dead) ;

INCLUDES *
sc_create_raster_nodata_includes(rt_raster r, int band) ;
void
sc_destroy_raster_nodata_includes(INCLUDES *dead) ;


EVALUATOR *
sc_create_raster_bands_evaluator(rt_raster raster, int *bands, int num_bands);
void
sc_destroy_raster_bands_evaluator(EVALUATOR *dead) ;


SPATIAL_COLLECTION *
sc_create_raster_wrapper(rt_raster raster, int owned, int *bands, int num_bands);
void
sc_destroy_raster_wrapper(SPATIAL_COLLECTION *dead);

SPATIAL_COLLECTION *
sc_create_raster_nodata_wrapper(rt_raster raster, int owned,
		int *bands, int num_bands, int nodata_band) ;
void
sc_destroy_raster_nodata_wrapper(SPATIAL_COLLECTION *dead) ;

SPATIAL_COLLECTION *
sc_create_raster_aligned_collection(SPATIAL_COLLECTION *wrapped, rt_raster alignTo);

void
sc_destroy_raster_aligned_collection(SPATIAL_COLLECTION *dead);


void
fit_raster_to_extent(GBOX *extent, rt_raster raster) ;

void
sync_raster_bands_to_value(VALUE *val, rt_raster empty_raster);

void
sc_sampling_engine(SPATIAL_COLLECTION *source,
		           rt_raster result,
		           VALUE *nodata_val) ;



#endif