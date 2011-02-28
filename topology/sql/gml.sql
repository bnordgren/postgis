-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
-- 
-- PostGIS - Spatial Types for PostgreSQL
-- http://postgis.refractions.net
--
-- Copyright (C) 2010, 2011 Sandro Santilli <strk@keybit.net>
--
-- This is free software; you can redistribute and/or modify it under
-- the terms of the GNU General Public Licence. See the COPYING file.
--
-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
--
-- Functions used for topology GML output
--
-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
--
-- Developed by Sandro Santilli <strk@keybit.net>
-- for Faunalia (http://www.faunalia.it) with funding from
-- Regione Toscana - Sistema Informativo per la Gestione del Territorio
-- e dell' Ambiente [RT-SIGTA].
-- For the project: "Sviluppo strumenti software per il trattamento di dati
-- geografici basati su QuantumGIS e Postgis (CIG 0494241492)"
--
-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

--{
--
-- INTERNAL FUNCTION
-- text _AsGMLNode(id, point, nsprefix, precision, options, idprefix, gmlver)
--
-- }{
CREATE OR REPLACE FUNCTION topology._AsGMLNode(id int, point geometry,
  nsprefix_in text, prec int, options int,
  idprefix text, gmlver int)
  RETURNS text
AS
$$
DECLARE
  nsprefix text;
  gml text;
BEGIN

  nsprefix := 'gml:';
  IF NOT nsprefix_in IS NULL THEN
    IF nsprefix_in = '' THEN
      nsprefix = nsprefix_in;
    ELSE
      nsprefix = nsprefix_in || ':';
    END IF;
  END IF;

  gml := '<' || nsprefix || 'Node ' || nsprefix
    || 'id="' || idprefix || 'N' || id || '"';
  IF point IS NOT NULL THEN
    gml = gml || '>'
              || '<' || nsprefix || 'pointProperty>'
              || ST_AsGML(gmlver, point, prec, options, nsprefix_in)
              || '</' || nsprefix || 'pointProperty>'
              || '</' || nsprefix || 'Node>';
  ELSE
    gml = gml || '/>';
  END IF;
  RETURN gml;
END
$$
LANGUAGE 'plpgsql';
--} _AsGMLNode(id, point, nsprefix, precision, options, idprefix, gmlVersion)

--{
--
-- INTERNAL FUNCTION
-- text _AsGMLEdge(edge_id, start_node, end_node, line, visitedTable,
--                 nsprefix, precision, options, idprefix, gmlVersion)
--
-- }{
CREATE OR REPLACE FUNCTION topology._AsGMLEdge(edge_id int, start_node int,
    end_node int, line geometry, visitedTable regclass, nsprefix_in text,
    prec int, options int, idprefix text, gmlver int)
  RETURNS text
AS
$$
DECLARE
  visited bool;
  nsprefix text;
  gml text;
BEGIN

  nsprefix := 'gml:';
  IF nsprefix_in IS NOT NULL THEN
    IF nsprefix_in = '' THEN
      nsprefix = nsprefix_in;
    ELSE
      nsprefix = nsprefix_in || ':';
    END IF;
  END IF;

  gml := '<' || nsprefix || 'Edge ' || nsprefix
    || 'id="' || idprefix || 'E' || edge_id || '">';

  -- Start node
  gml = gml || '<' || nsprefix || 'directedNode orientation="-"';
  -- Do visited bookkeeping if visitedTable was given
  visited = NULL;
  IF visitedTable IS NOT NULL THEN
    EXECUTE 'SELECT true FROM '
            || visitedTable::text
            || ' WHERE element_type = 1 AND element_id = '
            || start_node LIMIT 1 INTO visited;
    IF visited IS NOT NULL THEN
      gml = gml || ' xlink:href="#' || idprefix || 'N' || start_node || '" />';
    ELSE
      -- Mark as visited 
      EXECUTE 'INSERT INTO ' || visitedTable::text
        || '(element_type, element_id) VALUES (1, '
        || start_node || ')';
    END IF;
  END IF;
  IF visited IS NULL THEN
    gml = gml || '>';
    gml = gml || topology._AsGMLNode(start_node, NULL, nsprefix_in,
                                     prec, options, idprefix, gmlver);
    gml = gml || '</' || nsprefix || 'directedNode>';
  END IF;

  -- End node
  gml = gml || '<' || nsprefix || 'directedNode';
  -- Do visited bookkeeping if visitedTable was given
  visited = NULL;
  IF visitedTable IS NOT NULL THEN
    EXECUTE 'SELECT true FROM '
            || visitedTable::text
            || ' WHERE element_type = 1 AND element_id = '
            || end_node LIMIT 1 INTO visited;
    IF visited IS NOT NULL THEN
      gml = gml || ' xlink:href="#' || idprefix || 'N' || end_node || '" />';
    ELSE
      -- Mark as visited 
      EXECUTE 'INSERT INTO ' || visitedTable::text
        || '(element_type, element_id) VALUES (1, '
        || end_node || ')';
    END IF;
  END IF;
  IF visited IS NULL THEN
    gml = gml || '>';
    gml = gml || topology._AsGMLNode(end_node, NULL, nsprefix_in,
                                     prec, options, idprefix, gmlver);
    gml = gml || '</' || nsprefix || 'directedNode>';
  END IF;

  IF line IS NOT NULL THEN
    gml = gml || '<' || nsprefix || 'curveProperty>'
              || ST_AsGML(gmlver, line, prec, options, nsprefix_in)
              || '</' || nsprefix || 'curveProperty>';
  END IF;

  gml = gml || '</' || nsprefix || 'Edge>';

  RETURN gml;
END
$$
LANGUAGE 'plpgsql';
--} _AsGMLEdge(id, start_node, end_node, line, visitedTable, nsprefix, precision, options, idprefix, gmlver)

--{
--
-- INTERNAL FUNCTION
-- text _AsGMLFace(toponame, face_id, visitedTable,
--                 nsprefix, precision, options, idprefix, gmlVersion)
--
-- }{
CREATE OR REPLACE FUNCTION topology._AsGMLFace(toponame text, face_id int, 
    visitedTable regclass, nsprefix_in text,
    prec int, options int, idprefix text, gmlver int)
  RETURNS text
AS
$$
DECLARE
  visited bool;
  nsprefix text;
  gml text;
  rec RECORD;
  rec2 RECORD;
  bounds geometry;
  side int;
BEGIN

  nsprefix := 'gml:';
  IF nsprefix_in IS NOT NULL THEN
    IF nsprefix_in = '' THEN
      nsprefix = nsprefix_in;
    ELSE
      nsprefix = nsprefix_in || ':';
    END IF;
  END IF;

  gml := '<' || nsprefix || 'Face ' || nsprefix
    || 'id="' || idprefix || 'F' || face_id || '">';

  -- Construct the face geometry, then for each polygon:
  FOR rec IN SELECT (ST_DumpRings((ST_Dump(ST_ForceRHR(
    topology.ST_GetFaceGeometry(toponame, face_id)))).geom)).*
  LOOP

      -- Contents of a directed face are the list of edges
      -- that cover the specific ring
      bounds = ST_Boundary(rec.geom);

      FOR rec2 IN EXECUTE
        'SELECT e.*, ST_Line_Locate_Point('
        || quote_literal(bounds::text)
        || ', ST_Line_Interpolate_Point(e.geom, 0.2)) as pos FROM '
        || quote_ident(toponame)
        || '.edge e WHERE ( e.left_face = ' || face_id
        || ' OR e.right_face = ' || face_id
        || ') AND ST_Covers('
        || quote_literal(bounds::text)
        || ', e.geom) ORDER BY pos'
      LOOP

        gml = gml || '<' || nsprefix || 'directedEdge';

        -- if this edge goes in same direction to the
        --       ring bounds, make it with negative orientation
        SELECT DISTINCT (ST_Dump(
                          ST_SharedPaths(rec2.geom, bounds))
                        ).path[1] into side;
        IF side = 1 THEN -- edge goes in same direction
          gml = gml || ' orientation="-"';
        END IF;

        -- Do visited bookkeeping if visitedTable was given
        IF visitedTable IS NOT NULL THEN

          EXECUTE 'SELECT true FROM '
            || visitedTable::text
            || ' WHERE element_type = 2 AND element_id = '
            || rec2.edge_id LIMIT 1 INTO visited;
          IF visited THEN
            -- Use xlink:href if visited
            gml = gml || ' xlink:href="#' || idprefix || 'E'
                      || rec2.edge_id || '" />';
            CONTINUE;
          ELSE
            -- Mark as visited otherwise
            EXECUTE 'INSERT INTO ' || visitedTable::text
              || '(element_type, element_id) VALUES (2, '
              || rec2.edge_id || ')';
          END IF;

        END IF;

        gml = gml || '>';

        gml = gml || topology._AsGMLEdge(rec2.edge_id, rec2.start_node,
                                        rec2.end_node, rec2.geom,
                                        visitedTable, nsprefix_in,
                                        prec, options, idprefix, gmlver);
        gml = gml || '</' || nsprefix || 'directedEdge>';

      END LOOP;
    END LOOP;

  gml = gml || '</' || nsprefix || 'Face>';

  RETURN gml;
END
$$
LANGUAGE 'plpgsql';
--} _AsGMLFace(toponame, id, visitedTable, nsprefix, precision, options, idprefix, gmlver)

--{
--
-- API FUNCTION
--
-- text AsGML(TopoGeometry, nsprefix, precision, options, visitedTable, idprefix, gmlver)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry,
    nsprefix_in text, precision_in int, options_in int, visitedTable regclass,
    idprefix text, gmlver int)
  RETURNS text
AS
$$
DECLARE
  nsprefix text;
  precision int;
  options int;
  visited bool;
  toponame text;
  gml text;
  sql text;
  rec RECORD;
  rec2 RECORD;
  --bounds geometry;
  side int;
BEGIN

  nsprefix := 'gml:';
  IF nsprefix_in IS NOT NULL THEN
    IF nsprefix_in = '' THEN
      nsprefix = nsprefix_in;
    ELSE
      nsprefix = nsprefix_in || ':';
    END IF;
  END IF;

  precision := 15;
  IF precision_in IS NOT NULL THEN
    precision = precision_in;
  END IF;

  options := 1;
  IF options_in IS NOT NULL THEN
    options = options_in;
  END IF;

  -- Get topology name (for subsequent queries)
  SELECT name FROM topology.topology into toponame
              WHERE id = tg.topology_id;

  -- Puntual TopoGeometry
  IF tg.type = 1 THEN
    gml = '<' || nsprefix || 'TopoPoint>';
    -- For each defining node, print a directedNode
    FOR rec IN  EXECUTE 'SELECT r.element_id, n.geom from '
      || quote_ident(toponame) || '.relation r LEFT JOIN '
      || quote_ident(toponame) || '.node n ON (r.element_id = n.node_id)'
      || ' WHERE r.layer_id = ' || tg.layer_id
      || ' AND r.topogeo_id = ' || tg.id
    LOOP
      gml = gml || '<' || nsprefix || 'directedNode';
      -- Do visited bookkeeping if visitedTable was given
      IF visitedTable IS NOT NULL THEN
        EXECUTE 'SELECT true FROM '
                || visitedTable::text
                || ' WHERE element_type = 1 AND element_id = '
                || rec.element_id LIMIT 1 INTO visited;
        IF visited IS NOT NULL THEN
          gml = gml || ' xlink:href="#' || idprefix || 'N' || rec.element_id || '" />';
          CONTINUE;
        ELSE
          -- Mark as visited 
          EXECUTE 'INSERT INTO ' || visitedTable::text
            || '(element_type, element_id) VALUES (1, '
            || rec.element_id || ')';
        END IF;
      END IF;
      gml = gml || '>';
      gml = gml || topology._AsGMLNode(rec.element_id, rec.geom, nsprefix_in, precision, options, idprefix, gmlver);
      gml = gml || '</' || nsprefix || 'directedNode>';
    END LOOP;
    gml = gml || '</' || nsprefix || 'TopoPoint>';
    RETURN gml;

  ELSIF tg.type = 2 THEN -- lineal
    gml = '<' || nsprefix || 'TopoCurve>';

    FOR rec IN SELECT (ST_Dump(topology.Geometry(tg))).geom
    LOOP
      FOR rec2 IN EXECUTE
        'SELECT e.*, ST_Line_Locate_Point('
        || quote_literal(rec.geom::text)
        || ', ST_Line_Interpolate_Point(e.geom, 0.2)) as pos FROM '
        || quote_ident(toponame)
        || '.edge e WHERE ST_Covers('
        || quote_literal(rec.geom::text)
        || ', e.geom) ORDER BY pos'
        -- TODO: add relation to the conditional, to reduce load ?
      LOOP

        gml = gml || '<' || nsprefix || 'directedEdge';

        -- if this edge goes in opposite direction to the
        --       line, make it with negative orientation
        SELECT DISTINCT (ST_Dump(
                          ST_SharedPaths(rec2.geom, rec.geom))
                        ).path[1] into side;
        IF side = 2 THEN -- edge goes in opposite direction
          gml = gml || ' orientation="-"';
        END IF;

        -- Do visited bookkeeping if visitedTable was given
        IF visitedTable IS NOT NULL THEN

          EXECUTE 'SELECT true FROM '
            || visitedTable::text
            || ' WHERE element_type = 2 AND element_id = '
            || rec2.edge_id LIMIT 1 INTO visited;
          IF visited THEN
            -- Use xlink:href if visited
            gml = gml || ' xlink:href="#' || idprefix || 'E' || rec2.edge_id || '" />';
            CONTINUE;
          ELSE
            -- Mark as visited otherwise
            EXECUTE 'INSERT INTO ' || visitedTable::text
              || '(element_type, element_id) VALUES (2, '
              || rec2.edge_id || ')';
          END IF;

        END IF;


        gml = gml || '>';

        gml = gml || topology._AsGMLEdge(rec2.edge_id,
                                        rec2.start_node,
                                        rec2.end_node, rec2.geom,
                                        visitedTable,
                                        nsprefix_in, precision,
                                        options, idprefix, gmlver);


        gml = gml || '</' || nsprefix || 'directedEdge>';
      END LOOP;
    END LOOP;

    gml = gml || '</' || nsprefix || 'TopoCurve>';
    return gml;

  ELSIF tg.type = 3 THEN -- areal
    gml = '<' || nsprefix || 'TopoSurface>';

    -- For each defining face, print a directedFace
    FOR rec IN  EXECUTE 'SELECT f.face_id from '
      || quote_ident(toponame) || '.relation r LEFT JOIN '
      || quote_ident(toponame) || '.face f ON (r.element_id = f.face_id)'
      || ' WHERE r.layer_id = ' || tg.layer_id
      || ' AND r.topogeo_id = ' || tg.id
    LOOP
      gml = gml || '<' || nsprefix || 'directedFace';
      -- Do visited bookkeeping if visitedTable was given
      IF visitedTable IS NOT NULL THEN
        EXECUTE 'SELECT true FROM '
                || visitedTable::text
                || ' WHERE element_type = 3 AND element_id = '
                || rec.face_id LIMIT 1 INTO visited;
        IF visited IS NOT NULL THEN
          gml = gml || ' xlink:href="#' || idprefix || 'F' || rec.face_id || '" />';
          CONTINUE;
        ELSE
          -- Mark as visited 
          EXECUTE 'INSERT INTO ' || visitedTable::text
            || '(element_type, element_id) VALUES (3, '
            || rec.face_id || ')';
        END IF;
      END IF;
      gml = gml || '>';
      gml = gml || topology._AsGMLFace(toponame, rec.face_id, visitedTable,
                                       nsprefix_in, precision,
                                       options, idprefix, gmlver);
      gml = gml || '</' || nsprefix || 'directedFace>';
    END LOOP;
    gml = gml || '</' || nsprefix || 'TopoSurface>';
    RETURN gml;

  ELSIF tg.type = 4 THEN -- collection
    RAISE EXCEPTION 'Collection TopoGeometries are not supported by AsGML';

  END IF;
	

  RETURN gml;
	
END
$$
LANGUAGE 'plpgsql';
--} AsGML(TopoGeometry, nsprefix, precision, options, visitedTable, idprefix, gmlver)

--{
--
-- API FUNCTION
--
-- text AsGML(TopoGeometry, nsprefix, precision, options, visitedTable,
--            idprefix)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry,
    nsprefix text, prec int, options int, visitedTable regclass, idprefix text)
  RETURNS text
AS
$$
 SELECT topology.AsGML($1, $2, $3, $4, $5, $6, 3);
$$
LANGUAGE 'sql';
--} AsGML(TopoGeometry, nsprefix, precision, options, visitedTable, idprefix)

--{
--
-- API FUNCTION 
--
-- text AsGML(TopoGeometry, nsprefix, precision, options, visitedTable)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry,
    nsprefix text, prec int, options int, vis regclass)
  RETURNS text AS
$$
 SELECT topology.AsGML($1, $2, $3, $4, $5, '');
$$ LANGUAGE 'sql';
-- } AsGML(TopoGeometry, nsprefix, precision, options)


--{
--
-- API FUNCTION 
--
-- text AsGML(TopoGeometry, nsprefix, precision, options)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry,
    nsprefix text, prec int, opts int)
  RETURNS text AS
$$
 SELECT topology.AsGML($1, $2, $3, $4, NULL);
$$ LANGUAGE 'sql';
-- } AsGML(TopoGeometry, nsprefix, precision, options)

--{
--
-- API FUNCTION 
--
-- text AsGML(TopoGeometry, nsprefix)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry, nsprefix text)
  RETURNS text AS
$$
 SELECT topology.AsGML($1, $2, 15, 1, NULL);
$$ LANGUAGE 'sql';
-- } AsGML(TopoGeometry, nsprefix)

--{
--
-- API FUNCTION
--
-- text AsGML(TopoGeometry, visited_table)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry, visitedTable regclass)
  RETURNS text AS
$$
 SELECT topology.AsGML($1, 'gml', 15, 1, $2);
$$ LANGUAGE 'sql';
-- } AsGML(TopoGeometry, visited_table)

--{
--
-- API FUNCTION
--
-- text AsGML(TopoGeometry, visited_table, nsprefix)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry, visitedTable regclass, nsprefix text)
  RETURNS text AS
$$
 SELECT topology.AsGML($1, $3, 15, 1, $2);
$$ LANGUAGE 'sql';
-- } AsGML(TopoGeometry, visited_table, nsprefix)


--{
--
-- API FUNCTION
--
-- text AsGML(TopoGeometry)
--
-- }{
CREATE OR REPLACE FUNCTION topology.AsGML(tg topology.TopoGeometry)
  RETURNS text AS
$$
 SELECT topology.AsGML($1, 'gml');
$$ LANGUAGE 'sql';
-- } AsGML(TopoGeometry)

