/**********************************************************************
 * $Id:$
 *
 * PostGIS - Spatial Types for PostgreSQL
 * http://postgis.refractions.net
 * Copyright 2009 Oslandia
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU General Public Licence. See the COPYING file.
 *
 **********************************************************************/

/**
* @file KML input routines.
* Ability to parse KML geometry fragment and to return an LWGEOM
* or an error message.
*
* KML version supported: 2.2.0
* Cf: <http://www.opengeospatial.org/standards/kml>
*
* Known limitations related to 3D:
*  - Not support kml:Model geometries
*  - Don't handle kml:extrude attribute
*
* Written by Olivier Courtin - Oslandia
*
**********************************************************************/


#include "postgres.h"
#include "lwgeom_pg.h"
#include "liblwgeom.h"


#if HAVE_LIBXML2
#include <libxml/tree.h> 
#include <libxml/parser.h> 


/* 
TODO:
	- OGC:LonLat84_5773 explicit support (rather than EPSG:4326)
	- Don't return a GEOMETRYCOLLECTION if a MULTI one is enough
	- altitudeModeGroup relativeToGround Z Altitude
	  computation upon Geoid
*/


Datum geom_from_kml(PG_FUNCTION_ARGS);
static LWGEOM* parse_kml(xmlNodePtr xnode, bool *hasz);

#define KML_NS		((char *) "http://www.opengis.net/kml/2.2")


/**
 * Ability to parse KML geometry fragment and to return an LWGEOM
 * or an error message.
 */
PG_FUNCTION_INFO_V1(geom_from_kml);
Datum geom_from_kml(PG_FUNCTION_ARGS)
{
	PG_LWGEOM *geom, *geom2d;
	xmlDocPtr xmldoc;
	text *xml_input; 
	LWGEOM *lwgeom;
        int xml_size;
	uchar *srl;
	char *xml;
        size_t size=0;
	bool hasz=true;
	xmlNodePtr xmlroot=NULL;


	/* Get the KML stream */
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();
	xml_input = PG_GETARG_TEXT_P(0);

	xml_size = VARSIZE(xml_input) - VARHDRSZ; 	/* actual letters */
        xml = palloc(xml_size + 1); 			/* +1 for null */
	memcpy(xml, VARDATA(xml_input), xml_size);
	xml[xml_size] = 0; 				/* null term */

	/* Begin to Parse XML doc */
        xmlInitParser();
        xmldoc = xmlParseMemory(xml, xml_size);
        if (!xmldoc || (xmlroot = xmlDocGetRootElement(xmldoc)) == NULL) {
	        xmlFreeDoc(xmldoc);
	        xmlCleanupParser();
		lwerror("invalid KML representation");
	}

	lwgeom = parse_kml(xmlroot, &hasz);
	lwgeom->bbox = lwgeom_compute_box2d(lwgeom);
	geom = pglwgeom_serialize(lwgeom);
	lwgeom_release(lwgeom);

	xmlFreeDoc(xmldoc);
	xmlCleanupParser();

	/* KML geometries could be either 2 or 3D
	 *
	 * So we deal with 3D in all structures allocation, and flag hasz
	 * to false if we met once a missing Z dimension
	 * In this case, we force recursive 2D.
	 */
	if (!hasz) {
		srl = lwalloc(VARSIZE(geom));
        	lwgeom_force2d_recursive(SERIALIZED_FORM(geom), srl, &size);
        	geom2d = PG_LWGEOM_construct(srl, pglwgeom_getSRID(geom),
                                     lwgeom_hasBBOX(geom->type));
		lwfree(geom);
		geom = geom2d;
	}

	PG_RETURN_POINTER(geom);
}
			

/**
 * Return false if current element namespace is not a KML one
 * Return true otherwise.
 */
static bool is_kml_namespace(xmlNodePtr xnode, bool is_strict)
{
	xmlNsPtr *ns, *p;
	  
	ns = xmlGetNsList(xnode->doc, xnode);
	/*
	 * If no namespace is available we could return true anyway
	 * (because we work only on KML fragment, we don't want to 
	 *  'oblige' to add namespace on the geometry root node)
	 */
	if (ns == NULL) return !is_strict;

	for (p=ns ; *p ; p++) {
		if ((*p)->href == NULL) continue;
		if (!strcmp((char *) (*p)->href, KML_NS)) {
			if (	(*p)->prefix == NULL ||
				!xmlStrcmp(xnode->ns->prefix, (*p)->prefix)) {

				xmlFree(ns);
				return true;
			}
		}
	}

	xmlFree(ns);
	return false;
}


/**
 * Retrieve a KML propertie from a node or NULL otherwise
 * Respect namespaces if presents in the node element
 */
static xmlChar *kmlGetProp(xmlNodePtr xnode, xmlChar *prop)
{
	xmlChar *value;

	if (!is_kml_namespace(xnode, true))
		return xmlGetProp(xnode, prop);

	value = xmlGetNsProp(xnode, prop, (xmlChar *) KML_NS);

	/* In last case try without explicit namespace */
	if (value == NULL) value = xmlGetNoNsProp(xnode, prop);

	return value;
}


/**
 * Parse a string supposed to be a double
 */
static double parse_kml_double(char *d, bool space_before, bool space_after)
{
	char *p;
	int st;
	enum states {
		INIT     	= 0,
		NEED_DIG  	= 1,
		DIG	  	= 2,
		NEED_DIG_DEC 	= 3,
		DIG_DEC 	= 4,
		EXP	 	= 5,
		NEED_DIG_EXP 	= 6,
		DIG_EXP 	= 7,
		END 		= 8
	};

	/*
	 * Double pattern
	 * [-|\+]?[0-9]+(\.)?([0-9]+)?([Ee](\+|-)?[0-9]+)?
	 * We could also meet spaces before and/or after
	 * this pattern upon parameters
	 */

	if (space_before) while (isspace(*d)) d++;
	for (st = INIT, p = d ; *p ; p++) {

		if (isdigit(*p)) {
				if (st == INIT || st == NEED_DIG) 	st = DIG;
			else if (st == NEED_DIG_DEC) 			st = DIG_DEC;
			else if (st == NEED_DIG_EXP || st == EXP) 	st = DIG_EXP;
			else if (st == DIG || st == DIG_DEC || st == DIG_EXP);
			else lwerror("invalid KML representation"); 
		} else if (*p == '.') {
			if      (st == DIG) 				st = NEED_DIG_DEC;
			else    lwerror("invalid KML representation"); 
		} else if (*p == '-' || *p == '+') {
			if      (st == INIT) 				st = NEED_DIG;
			else if (st == EXP) 				st = NEED_DIG_EXP;
			else    lwerror("invalid KML representation"); 
		} else if (*p == 'e' || *p == 'E') {
			if      (st == DIG || st == DIG_DEC) 		st = EXP;
			else    lwerror("invalid KML representation"); 
		} else if (isspace(*p)) {
			if (!space_after) lwerror("invalid KML representation");  
			if (st == DIG || st == DIG_DEC || st == DIG_EXP)st = END;
			else if (st == NEED_DIG_DEC)			st = END;
			else if (st == END);
			else    lwerror("invalid KML representation");
		} else  lwerror("invalid KML representation");
     	}	       

	if (st != DIG && st != NEED_DIG_DEC && st != DIG_DEC && st != DIG_EXP && st != END)
		lwerror("invalid KML representation");

	return atof(d);
}


/**
 * Parse kml:coordinates
 */
static POINTARRAY* parse_kml_coordinates(xmlNodePtr xnode, bool *hasz)
{
        xmlChar *kml_coord;
        bool digit, found;
        DYNPTARRAY *dpa;
        POINTARRAY *pa;
        int kml_dims;
        char *p, *q;
        POINT4D pt;
        uchar dims=0;

	if (xnode == NULL) lwerror("invalid KML representation");
	
	for (found = false ; xnode != NULL ; xnode = xnode->next) {
		if (xnode->type != XML_ELEMENT_NODE) continue;
		if (!is_kml_namespace(xnode, false)) continue;
		if (strcmp((char *) xnode->name, "coordinates")) continue;

		found = true;
		break;
	}
	if (!found) lwerror("invalid KML representation");

        /* We begin to retrieve coordinates string */
        kml_coord = xmlNodeGetContent(xnode);
        p = (char *) kml_coord;

        /* KML coordinates pattern:     x1,y1 x2,y2 
         *                              x1,y1,z1 x2,y2,z2
	 */

        /* Now we create PointArray from coordinates values */
        TYPE_SETZM(dims, 1, 0);
        dpa = dynptarray_create(1, dims);

        for (q = p, kml_dims=0, digit = false ; *p ; p++) {

                if (isdigit(*p)) digit = true;  /* One state parser */

                /* Coordinate Separator */
                if (*p == ',') {
                        *p = '\0';
                        kml_dims++;

                        if (*(p+1) == '\0') lwerror("invalid KML representation");

                        if      (kml_dims == 1) pt.x = parse_kml_double(q, true, true);
                        else if (kml_dims == 2) pt.y = parse_kml_double(q, true, true);
                        q = p+1;

                /* Tuple Separator (or end string) */
                } else if (digit && (isspace(*p) || *(p+1) == '\0')) {
                        if (isspace(*p)) *p = '\0';
                        kml_dims++;

                        if (kml_dims < 2 || kml_dims > 3)
                                lwerror("invalid KML representation");

                        if (kml_dims == 3)
                                pt.z = parse_kml_double(q, true, true);
                        else {
                                pt.y = parse_kml_double(q, true, true);
                                *hasz = false;
                        }

                        dynptarray_addPoint4d(dpa, &pt, 0);
                        digit = false;
                        q = p+1;
                        kml_dims = 0;

                }
        }

        xmlFree(kml_coord);
        pa = ptarray_clone(dpa->pa);
        lwfree(dpa);

        return pa;
}


/**
 * Parse KML point
 */
static LWGEOM* parse_kml_point(xmlNodePtr xnode, bool *hasz)
{
	POINTARRAY *pa;

	if (xnode->children == NULL) lwerror("invalid KML representation");
	pa = parse_kml_coordinates(xnode->children, hasz);
	if (pa->npoints != 1) lwerror("invalid KML representation");

	return (LWGEOM *) lwpoint_construct(4326, NULL, pa);
}


/**
 * Parse KML lineString 
 */
static LWGEOM* parse_kml_line(xmlNodePtr xnode, bool *hasz)
{
	POINTARRAY *pa;

	if (xnode->children == NULL) lwerror("invalid KML representation");
	pa = parse_kml_coordinates(xnode->children, hasz);
	if (pa->npoints < 2) lwerror("invalid KML representation");

	return (LWGEOM *) lwline_construct(4326, NULL, pa);
}


/**
 * Parse KML Polygon
 */
static LWGEOM* parse_kml_polygon(xmlNodePtr xnode, bool *hasz)
{
	int ring;
	xmlNodePtr xa, xb;
	POINTARRAY **ppa = NULL;

	for (xa = xnode->children ; xa != NULL ; xa = xa->next) {

		/* Polygon/outerBoundaryIs */
		if (xa->type != XML_ELEMENT_NODE) continue;
		if (!is_kml_namespace(xa, false)) continue;
		if (strcmp((char *) xa->name, "outerBoundaryIs")) continue;
	       
		for (xb = xa->children ; xb != NULL ; xb = xb->next) {

			if (xb->type != XML_ELEMENT_NODE) continue;
			if (!is_kml_namespace(xb, false)) continue;
			if (strcmp((char *) xb->name, "LinearRing")) continue;

			ppa = (POINTARRAY**) lwalloc(sizeof(POINTARRAY*));
			ppa[0] = parse_kml_coordinates(xb->children, hasz);

			if (ppa[0]->npoints < 4
				|| (!*hasz && !ptarray_isclosed2d(ppa[0]))
				||  (*hasz && !ptarray_isclosed3d(ppa[0])))
				lwerror("invalid KML representation");
		}
	}

	for (ring=1, xa = xnode->children ; xa != NULL ; xa = xa->next) {

		/* Polygon/innerBoundaryIs */
		if (xa->type != XML_ELEMENT_NODE) continue;
		if (!is_kml_namespace(xa, false)) continue;
		if (strcmp((char *) xa->name, "innerBoundaryIs")) continue;
		
		for (xb = xa->children ; xb != NULL ; xb = xb->next) {

			if (xb->type != XML_ELEMENT_NODE) continue;
			if (!is_kml_namespace(xb, false)) continue;
			if (strcmp((char *) xb->name, "LinearRing")) continue;

			ppa = (POINTARRAY**) lwrealloc((POINTARRAY *) ppa,
				sizeof(POINTARRAY*) * (ring + 1));
			ppa[ring] = parse_kml_coordinates(xb->children, hasz);

			if (ppa[ring]->npoints < 4
				|| (!*hasz && !ptarray_isclosed2d(ppa[ring]))
				||  (*hasz && !ptarray_isclosed3d(ppa[ring])))
				lwerror("invalid KML representation");

			ring++;
		}
	}
			
	/* Exterior Ring is mandatory */
	if (ppa == NULL || ppa[0] == NULL) lwerror("invalid KML representation");

	return (LWGEOM *) lwpoly_construct(4326, NULL, ring, ppa);
}


/**
 * Parse KML MultiGeometry
 */
static LWGEOM* parse_kml_multi(xmlNodePtr xnode, bool *hasz) 
{
	LWGEOM *geom;
	xmlNodePtr xa;

 	geom = (LWGEOM *)lwcollection_construct_empty(4326, 1, 0);

	for (xa = xnode->children ; xa != NULL ; xa = xa->next) {

		if (xa->type != XML_ELEMENT_NODE) continue;
		if (!is_kml_namespace(xa, false)) continue;

		if (	   !strcmp((char *) xa->name, "Point")
		 	|| !strcmp((char *) xa->name, "LineString")
		 	|| !strcmp((char *) xa->name, "Polygon")
		 	|| !strcmp((char *) xa->name, "MultiGeometry")) {

			if (xa->children == NULL) break;
			geom = lwcollection_add((LWCOLLECTION *)geom, -1,
				 parse_kml(xa, hasz));
		}
	}

	return geom;
}


/**
 * Parse KML 
 */
static LWGEOM* parse_kml(xmlNodePtr xnode, bool *hasz) 
{
	xmlNodePtr xa = xnode;

	while (xa != NULL && (xa->type != XML_ELEMENT_NODE
			|| !is_kml_namespace(xa, false))) xa = xa->next;

	if (xa == NULL) lwerror("invalid KML representation");

	if (!strcmp((char *) xa->name, "Point"))
		return parse_kml_point(xa, hasz);

	if (!strcmp((char *) xa->name, "LineString"))
		return parse_kml_line(xa, hasz);

	if (!strcmp((char *) xa->name, "Polygon"))
		return parse_kml_polygon(xa, hasz);

	if (!strcmp((char *) xa->name, "MultiGeometry"))
		return parse_kml_multi(xa, hasz);
	
	lwerror("invalid KML representation");
	return NULL; /* Never reach */
}

#endif /* if HAVE_LIBXML2 */