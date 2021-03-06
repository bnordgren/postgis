/**********************************************************************
 * $Id: cu_print.c 6160 2010-11-01 01:28:12Z pramsey $
 *
 * PostGIS - Spatial Types for PostgreSQL
 * http://postgis.refractions.net
 *
 * Copyright (C) 2011 Sandro Santilli <strk@keybit.net>
 * Copyright (C) 2008 Paul Ramsey
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU General Public Licence. See the COPYING file.
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CUnit/Basic.h"
#include "CUnit/CUnit.h"

#include "liblwgeom_internal.h"
#include "cu_tester.h"


static LWGEOM* lwgeom_from_text(const char *str)
{
	LWGEOM_PARSER_RESULT r;
	if( LW_FAILURE == lwgeom_parse_wkt(&r, (char*)str, LW_PARSER_CHECK_NONE) )
		return NULL;
	return r.geom;
}

static char* lwgeom_to_text(const LWGEOM *geom)
{
	return lwgeom_to_wkt(geom, WKT_ISO, 8, NULL);
}

static void test_ptarray_append_point(void)
{
	LWLINE *line;
	char *wkt;
	POINT4D p;

	line = lwgeom_as_lwline(lwgeom_from_text("LINESTRING(0 0,1 1)"));
	p.x = 1; 
	p.y = 1;
	ptarray_append_point(line->points, &p, LW_TRUE);
	wkt = lwgeom_to_text(lwline_as_lwgeom(line));
	CU_ASSERT_STRING_EQUAL(wkt,"LINESTRING(0 0,1 1,1 1)");
	lwfree(wkt);

	ptarray_append_point(line->points, &p, LW_FALSE);
	wkt = lwgeom_to_text(lwline_as_lwgeom(line));
	CU_ASSERT_STRING_EQUAL(wkt,"LINESTRING(0 0,1 1,1 1)");
	lwfree(wkt);

	lwline_free(line);
}

static void test_ptarray_append_ptarray(void)
{
	/* TODO: add tests here ! */
}

static void test_ptarray_locate_point(void)
{
        LWLINE *line;
	double loc, dist;
        POINT2D p;

        line = lwgeom_as_lwline(lwgeom_from_text("LINESTRING(0 3,20 4)"));

	p = getPoint2d(line->points, 0);
	loc = ptarray_locate_point(line->points, &p, &dist);
	CU_ASSERT_EQUAL(loc, 0);
	CU_ASSERT_EQUAL(dist, 0.0);

	p = getPoint2d(line->points, 1);
	loc = ptarray_locate_point(line->points, &p, &dist);
	CU_ASSERT_EQUAL(loc, 1);
	CU_ASSERT_EQUAL(dist, 0.0);

	p.x = 21; p.y = 4;
	loc = ptarray_locate_point(line->points, &p, &dist);
	CU_ASSERT_EQUAL(loc, 1);
	CU_ASSERT_EQUAL(dist, 1.0);

	p.x = 0; p.y = 2;
	loc = ptarray_locate_point(line->points, &p, &dist);
	CU_ASSERT_EQUAL(loc, 0);
	CU_ASSERT_EQUAL(dist, 1.0);

	lwline_free(line);
	line = lwgeom_as_lwline(lwgeom_from_text("LINESTRING(0 0,20 0,40 0)"));

	p.x = 20; p.y = 0;
	loc = ptarray_locate_point(line->points, &p, &dist);
	CU_ASSERT_EQUAL(loc, 0.5);
	CU_ASSERT_EQUAL(dist, 0.0);

	lwline_free(line);
	line = lwgeom_as_lwline(lwgeom_from_text("LINESTRING(-40 0,0 0,20 0,40 0)"));

	p.x = 20; p.y = 0;
	loc = ptarray_locate_point(line->points, &p, &dist);
	CU_ASSERT_EQUAL(loc, 0.75);
	CU_ASSERT_EQUAL(dist, 0.0);

	lwline_free(line);
}

static void test_ptarray_isccw(void)
{
        LWLINE *line;
	LWPOLY* poly;
	int ccw;

	/* clockwise rectangle */
        line = lwgeom_as_lwline(lwgeom_from_text("LINESTRING(0 0,0 10,10 10,10 0, 0 0)"));
	ccw = ptarray_isccw(line->points);
	CU_ASSERT_EQUAL(ccw, 0);
	lwline_free(line);

	/* clockwise triangle */
        line = lwgeom_as_lwline(lwgeom_from_text("LINESTRING(0 3,20 4,20 3, 0 3)"));
	ccw = ptarray_isccw(line->points);
	CU_ASSERT_EQUAL(ccw, 0);
	lwline_free(line);

	/* counterclockwise triangle */
        line = lwgeom_as_lwline(lwgeom_from_text("LINESTRING(0 3,20 3,20 4, 0 3)"));
	ccw = ptarray_isccw(line->points);
	CU_ASSERT_EQUAL(ccw, 1);
	lwline_free(line);

	/* counterclockwise narrow ring (see ticket #1302) */
        line = lwgeom_as_lwline(lwgeom_from_hexwkb("01020000000500000000917E9BA468294100917E9B8AEA284137894120A4682941C976BE9F8AEA2841B39ABE1FA46829415ACCC29F8AEA2841C976BE1FA4682941C976BE9F8AEA284100917E9BA468294100917E9B8AEA2841", LW_PARSER_CHECK_NONE));
	ccw = ptarray_isccw(line->points);
	CU_ASSERT_EQUAL(ccw, 1);
	lwline_free(line);

	/* clockwise narrow ring (see ticket #1302) */
	line = lwgeom_as_lwline(lwgeom_from_hexwkb("01020000000500000000917E9BA468294100917E9B8AEA2841C976BE1FA4682941C976BE9F8AEA2841B39ABE1FA46829415ACCC29F8AEA284137894120A4682941C976BE9F8AEA284100917E9BA468294100917E9B8AEA2841", LW_PARSER_CHECK_NONE));
	ccw = ptarray_isccw(line->points);
	CU_ASSERT_EQUAL(ccw, 0);
	lwline_free(line);

	/* Clockwise narrow ring (see ticket #1302) */
	poly = lwgeom_as_lwpoly(lwgeom_from_hexwkb("0103000000010000000500000000917E9BA468294100917E9B8AEA2841C976BE1FA4682941C976BE9F8AEA2841B39ABE1FA46829415ACCC29F8AEA284137894120A4682941C976BE9F8AEA284100917E9BA468294100917E9B8AEA2841", LW_PARSER_CHECK_NONE));
	ccw = ptarray_isccw(poly->rings[0]);
	CU_ASSERT_EQUAL(ccw, 0);
	lwpoly_free(poly);
}


/*
** Used by the test harness to register the tests in this file.
*/
CU_TestInfo ptarray_tests[] =
{
	PG_TEST(test_ptarray_append_point),
	PG_TEST(test_ptarray_append_ptarray),
	PG_TEST(test_ptarray_locate_point),
	PG_TEST(test_ptarray_isccw),
	CU_TEST_INFO_NULL
};
CU_SuiteInfo ptarray_suite = {"ptarray", NULL, NULL, ptarray_tests };
