#ifndef RT_COLLECTION_H
#define RT_COLLECTION_H

#include "spatial_collection.h"
#include "lwgeom_pg.h"
#include "rt_pg.h"
#include "fmgr.h"

SPATIAL_COLLECTION *
sc_create_pgraster_wrapper_nodata(rt_pgraster *pg_raster,
		int *bands, int num_bands, int nodata_band) ;

SPATIAL_COLLECTION *
sc_create_pgraster_wrapper(rt_pgraster *pg_raster, int *bands, int num_bands) ;

SPATIAL_COLLECTION *
sc_create_pglwgeom_wrapper(PG_LWGEOM *pg_geom, double inside, double outside) ;

int
getarg_bandlist(FunctionCallInfo fcinfo, int argnum,
		        rt_pgraster *raster, int **bands, int *num_bands) ;

rt_raster
rt_raster_new_inbox(GBOX *extent, rt_pgraster *grid_defn) ;



#endif
