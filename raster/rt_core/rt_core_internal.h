#include "rt_api.h"

#ifndef RT_CORE_INTERNAL_H
#define RT_CORE_INTERNAL_H

struct rt_extband_t {
    uint8_t bandNum;
    char* path; /* externally owned ? */
};

struct rt_band_t {
    rt_pixtype pixtype;
    int32_t offline;
    uint16_t width;
    uint16_t height;
    int32_t hasnodata; /* a flag indicating if this band contains nodata values */
    int32_t isnodata;   /* a flag indicating if this band is filled only with
                           nodata values */
    double nodataval; /* int will be converted ... */
    int32_t ownsData; /* XXX mloskot: its behaviour needs to be documented */

    union {
        void* mem; /* actual data, externally owned */
        struct rt_extband_t offline;
    } data;

};

struct rt_bandstats_t {
	double sample;
	uint32_t count;

	double min;
	double max;
	double sum;
	double mean;
	double stddev;

	double *values;
	int sorted; /* flag indicating that values is sorted ascending by value */
};

struct rt_histogram_t {
	uint32_t count;
	double percent;

	double min;
	double max;

	int inc_min;
	int inc_max;
};

struct rt_quantile_t {
	double quantile;
	double value;
};

struct rt_valuecount_t {
	double value;
	uint32_t count;
	double percent;
};

struct rt_reclassexpr_t {
	struct rt_reclassrange {
		double min;
		double max;
		int inc_min; /* include min */
		int inc_max; /* include max */
		int exc_min; /* exceed min */
		int exc_max; /* exceed max */
	} src, dst;
};

struct rt_raster_serialized_t {
    /*---[ 8 byte boundary ]---{ */
    uint32_t size; /* required by postgresql: 4 bytes */
    uint16_t version; /* format version (this is version 0): 2 bytes */
    uint16_t numBands; /* Number of bands: 2 bytes */

    /* }---[ 8 byte boundary ]---{ */
    double scaleX; /* pixel width: 8 bytes */

    /* }---[ 8 byte boundary ]---{ */
    double scaleY; /* pixel height: 8 bytes */

    /* }---[ 8 byte boundary ]---{ */
    double ipX; /* insertion point X: 8 bytes */

    /* }---[ 8 byte boundary ]---{ */
    double ipY; /* insertion point Y: 8 bytes */

    /* }---[ 8 byte boundary ]---{ */
    double skewX; /* skew about the X axis: 8 bytes */

    /* }---[ 8 byte boundary ]---{ */
    double skewY; /* skew about the Y axis: 8 bytes */

    /* }---[ 8 byte boundary ]--- */
    int32_t srid; /* Spatial reference id: 4 bytes */
    uint16_t width; /* pixel columns: 2 bytes */
    uint16_t height; /* pixel rows: 2 bytes */
};

/* NOTE: the initial part of this structure matches the layout
 *       of data in the serialized form version 0, starting
 *       from the numBands element
 */
struct rt_raster_t {
    uint32_t size;
    uint16_t version;

    /* Number of bands, all share the same dimension
     * and georeference */
    uint16_t numBands;

    /* Georeference (in projection units) */
    double scaleX; /* pixel width */
    double scaleY; /* pixel height */
    double ipX; /* geo x ordinate of the corner of upper-left pixel */
    double ipY; /* geo y ordinate of the corner of bottom-right pixel */
    double skewX; /* skew about the X axis*/
    double skewY; /* skew about the Y axis */

    int32_t srid; /* spatial reference id */
    uint16_t width; /* pixel columns - max 65535 */
    uint16_t height; /* pixel rows - max 65535 */
    rt_band *bands; /* actual bands */

};

/* WKT string representing each polygon in WKT format acompagned by its
correspoding value */
struct rt_geomval_t {
    int srid;
    double val;
    char * geom;
};

struct rt_gdaldriver_t {
    int idx;
    char *short_name;
    char *long_name;
		char *create_options;
};

#endif /* RT_CORE_INTERNAL_H */
