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
-- Functions used to populate a topology
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
-- AddNode(atopology, point)
--
-- Add a node primitive to a topology and get it's identifier.
-- Returns an existing node at the same location, if any.
--
-- When adding a _new_ node it checks for the existance of any
-- edge crossing the given point, raising an exception if found.
--
-- The newly added nodes have no containing face.
--
-- 
CREATE OR REPLACE FUNCTION topology.AddNode(varchar, geometry)
	RETURNS int
AS
$$
DECLARE
	atopology ALIAS FOR $1;
	apoint ALIAS FOR $2;
	nodeid int;
	rec RECORD;
BEGIN
	--
	-- Atopology and apoint are required
	-- 
	IF atopology IS NULL OR apoint IS NULL THEN
		RAISE EXCEPTION 'Invalid null argument';
	END IF;

	--
	-- Apoint must be a point
	--
	IF substring(geometrytype(apoint), 1, 5) != 'POINT'
	THEN
		RAISE EXCEPTION 'Node geometry must be a point';
	END IF;

	--
	-- Check if a coincident node already exists
	-- 
	-- We use index AND x/y equality
	--
	FOR rec IN EXECUTE 'SELECT node_id FROM '
		|| quote_ident(atopology) || '.node ' ||
		'WHERE geom && ' || quote_literal(apoint::text) || '::geometry'
		||' AND ST_X(geom) = ST_X('||quote_literal(apoint::text)||'::geometry)'
		||' AND ST_Y(geom) = ST_Y('||quote_literal(apoint::text)||'::geometry)'
	LOOP
		RETURN  rec.node_id;
	END LOOP;

	--
	-- Check if any edge crosses this node
	-- (endpoints are fine)
	--
	FOR rec IN EXECUTE 'SELECT edge_id FROM '
		|| quote_ident(atopology) || '.edge ' 
		|| 'WHERE ST_DWithin('
		|| quote_literal(apoint::text) 
		|| ', geom, 0) AND NOT ST_Equals('
		|| quote_literal(apoint::text)
		|| ', ST_StartPoint(geom)) AND NOT ST_Equals('
		|| quote_literal(apoint::text)
		|| ', ST_EndPoint(geom))'
	LOOP
		RAISE EXCEPTION 'An edge crosses the given node.';
	END LOOP;

	--
	-- Get new node id from sequence
	--
	FOR rec IN EXECUTE 'SELECT nextval(' ||
		quote_literal(
			quote_ident(atopology) || '.node_node_id_seq'
		) || ')'
	LOOP
		nodeid = rec.nextval;
	END LOOP;

	--
	-- Insert the new row
	--
	EXECUTE 'INSERT INTO ' || quote_ident(atopology)
		|| '.node(node_id, geom) 
		VALUES('||nodeid||','||quote_literal(apoint::text)||
		')';

	RETURN nodeid;
	
END
$$
LANGUAGE 'plpgsql' VOLATILE;
--} AddNode

--{
--
-- AddEdge(atopology, line)
--
-- Add an edge primitive to a topology and get it's identifier.
-- Edge endpoints will be added as nodes if missing.
-- Returns an existing edge at the same location, if any.
--
-- An exception is raised if the given line crosses an existing
-- node or interects with an existing edge on anything but endnodes.
--
-- The newly added edge has "universe" face on both sides
-- and links to itself as per next left/right edge.
-- Calling code is expected to do further linking.
--
-- 
CREATE OR REPLACE FUNCTION topology.AddEdge(varchar, geometry)
	RETURNS int
AS
$$
DECLARE
	atopology ALIAS FOR $1;
	aline ALIAS FOR $2;
	edgeid int;
	rec RECORD;
  ix geometry; 
BEGIN
	--
	-- Atopology and apoint are required
	-- 
	IF atopology IS NULL OR aline IS NULL THEN
		RAISE EXCEPTION 'Invalid null argument';
	END IF;

	--
	-- Aline must be a linestring
	--
	IF substring(geometrytype(aline), 1, 4) != 'LINE'
	THEN
		RAISE EXCEPTION 'Edge geometry must be a linestring';
	END IF;

	--
	-- Check there's no face registered in the topology
	--
	FOR rec IN EXECUTE 'SELECT count(face_id) FROM '
		|| quote_ident(atopology) || '.face '
		|| ' WHERE face_id != 0 LIMIT 1'
	LOOP
		IF rec.count > 0 THEN
			RAISE EXCEPTION 'AddEdge can only be used against topologies with no faces defined';
		END IF;
	END LOOP;

	--
	-- Check if the edge crosses an existing node
	--
	FOR rec IN EXECUTE 'SELECT node_id FROM '
		|| quote_ident(atopology) || '.node '
		|| 'WHERE ST_Crosses('
		|| quote_literal(aline::text) || '::geometry, geom'
		|| ')'
	LOOP
		RAISE EXCEPTION 'Edge crosses node %', rec.node_id;
	END LOOP;

	--
	-- Check if the edge intersects an existing edge
	-- on anything but endpoints
	--
	-- Following DE-9 Intersection Matrix represent
	-- the only relation we accept. 
	--
	--    F F 1
	--    F * *
	--    1 * 2
	--
	-- Example1: linestrings touching at one endpoint
	--    FF1 F00 102
	--    FF1 F** 1*2 <-- our match
	--
	-- Example2: linestrings touching at both endpoints
	--    FF1 F0F 1F2
	--    FF1 F** 1*2 <-- our match
	--
	FOR rec IN EXECUTE 'SELECT edge_id, geom, ST_Relate('
		|| quote_literal(aline::text)
		|| '::geometry, geom, 2) as im'
		|| ' FROM '
		|| quote_ident(atopology) || '.edge '
		|| 'WHERE '
		|| quote_literal(aline::text) || '::geometry && geom'

	LOOP

	  IF ST_RelateMatch(rec.im, 'FF1F**1*2') THEN
	    CONTINUE; -- no interior intersection
	  END IF;

	  -- Reuse an EQUAL edge (be it closed or not)
	  IF ST_RelateMatch(rec.im, '1FFF*FFF2') THEN
	      RAISE DEBUG 'Edge already known as %', rec.edge_id;
	      RETURN rec.edge_id;
	  END IF;

	  -- WARNING: the constructive operation might throw an exception
	  BEGIN
	    ix = ST_Intersection(rec.geom, aline);
	  EXCEPTION
	  WHEN OTHERS THEN
	    RAISE NOTICE 'Could not compute intersection between input edge (%) and edge % (%)', aline::text, rec.edge_id, rec.geom::text;
	  END;

	  RAISE EXCEPTION 'Edge intersects (not on endpoints) with existing edge % at or near point %', rec.edge_id, ST_AsText(ST_PointOnSurface(ix));

	END LOOP;

	--
	-- Get new edge id from sequence
	--
	FOR rec IN EXECUTE 'SELECT nextval(' ||
		quote_literal(
			quote_ident(atopology) || '.edge_data_edge_id_seq'
		) || ')'
	LOOP
		edgeid = rec.nextval;
	END LOOP;

	--
	-- Insert the new row
	--
	EXECUTE 'INSERT INTO '
		|| quote_ident(atopology)
		|| '.edge(edge_id, start_node, end_node, '
		|| 'next_left_edge, next_right_edge, '
		|| 'left_face, right_face, '
		|| 'geom) '
		|| ' VALUES('

		-- edge_id
		|| edgeid ||','

		-- start_node
		|| 'topology.addNode('
		|| quote_literal(atopology)
		|| ', ST_StartPoint('
		|| quote_literal(aline::text)
		|| ')) ,'

		-- end_node
		|| 'topology.addNode('
		|| quote_literal(atopology)
		|| ', ST_EndPoint('
		|| quote_literal(aline::text)
		|| ')) ,'

		-- next_left_edge
		|| -edgeid ||','

		-- next_right_edge
		|| edgeid ||','

		-- left_face
		|| '0,'

		-- right_face
		|| '0,'

		-- geom
		||quote_literal(aline::text)
		|| ')';

	RETURN edgeid;
	
END
$$
LANGUAGE 'plpgsql' VOLATILE;
--} AddEdge

--{
--
-- AddFace(atopology, poly, [<force_new>=true])
--
-- Add a face primitive to a topology and get it's identifier.
-- Returns an existing face at the same location, if any, unless
-- true is passed as the force_new argument
--
-- For a newly added face, its edges will be appropriately
-- linked (marked as left-face or right-face), and any contained
-- edges and nodes would also be marked as such.
--
-- When forcing re-registration of an existing face, no action will be
-- taken to deal with the face being substituted. Which means
-- a record about the old face and any record in the relation table
-- referencing the existing face will remain untouched, effectively
-- leaving the topology in a possibly invalid state.
-- It is up to the caller to deal with that.
--
-- The target topology is assumed to be valid (containing no
-- self-intersecting edges).
--
-- An exception is raised if:
--  o The polygon boundary is not fully defined by existing edges.
--  o The polygon overlaps an existing face.
--
-- 
CREATE OR REPLACE FUNCTION topology.AddFace(atopology varchar, apoly geometry, force_new boolean DEFAULT FALSE)
	RETURNS int
AS
$$
DECLARE
  bounds geometry;
  symdif geometry;
  faceid int;
  rec RECORD;
  rrec RECORD;
  relate text;
  right_edges int[];
  left_edges int[];
  all_edges geometry;
  old_faceid int;
  old_edgeid int;
  sql text;
  right_side bool;
  edgeseg geometry;
  p1 geometry;
  p2 geometry;
  p3 geometry;
  loc float8;
  segnum int;
  numsegs int;
BEGIN
  --
  -- Atopology and apoly are required
  -- 
  IF atopology IS NULL OR apoly IS NULL THEN
    RAISE EXCEPTION 'Invalid null argument';
  END IF;

  --
  -- Aline must be a polygon
  --
  IF substring(geometrytype(apoly), 1, 4) != 'POLY'
  THEN
    RAISE EXCEPTION 'Face geometry must be a polygon';
  END IF;

  for rrec IN SELECT (ST_DumpRings(ST_ForceRHR(apoly))).*
  LOOP -- {
    --
    -- Find all bounds edges, forcing right-hand-rule
    -- to know what's left and what's right...
    --
    bounds = ST_Boundary(rrec.geom);

    sql := 'SELECT e.geom, e.edge_id, '
      || 'e.left_face, e.right_face FROM '
      || quote_ident(atopology) || '.edge e, (SELECT '
      || quote_literal(bounds::text)
      || '::geometry as geom) r WHERE '
      || 'r.geom && e.geom'
    ;
    -- RAISE DEBUG 'SQL: %', sql;
    FOR rec IN EXECUTE sql
    LOOP -- {
      --RAISE DEBUG 'Edge % has bounding box intersection', rec.edge_id;

      -- Find first non-empty segment of the edge
      numsegs = ST_NumPoints(rec.geom);
      segnum = 1;
      WHILE segnum < numsegs LOOP
        p1 = ST_PointN(rec.geom, segnum);
        p2 = ST_PointN(rec.geom, segnum+1);
        IF ST_Distance(p1, p2) > 0 THEN
          EXIT;
        END IF;
        segnum = segnum + 1;
      END LOOP;

      IF segnum = numsegs THEN
        RAISE WARNING 'Edge % is collapsed', rec.edge_id;
        CONTINUE; -- we don't want to spend time on it
      END IF;

      edgeseg = ST_MakeLine(p1, p2);

      -- Skip non-covered edges
      IF NOT ST_Equals(p2, ST_EndPoint(rec.geom)) THEN
        IF NOT ( _ST_Intersects(bounds, p1) AND _ST_Intersects(bounds, p2) )
        THEN
          --RAISE DEBUG 'Edge % has points % and % not intersecting with ring bounds', rec.edge_id, st_astext(p1), st_astext(p2);
          CONTINUE;
        END IF;
      ELSE
        -- must be a 2-points only edge, let's use Covers (more expensive)
        IF NOT _ST_Covers(bounds, edgeseg) THEN
          --RAISE DEBUG 'Edge % is not covered by ring', rec.edge_id;
          CONTINUE;
        END IF;
      END IF;

      p3 = ST_StartPoint(bounds);
      IF ST_DWithin(edgeseg, p3, 0) THEN
        -- Edge segment covers ring endpoint, See bug #874
        loc = ST_Line_Locate_Point(edgeseg, p3);
        -- WARNING: this is as robust as length of edgeseg allows...
        IF loc > 0.9 THEN
          -- shift last point down 
          p2 = ST_Line_Interpolate_Point(edgeseg, loc - 0.1);
        ELSIF loc < 0.1 THEN
          -- shift first point up
          p1 = ST_Line_Interpolate_Point(edgeseg, loc + 0.1); 
        ELSE
          -- when ring start point is in between, we swap the points
          p3 = p1; p1 = p2; p2 = p3;
        END IF;
      END IF;

      right_side = ST_Line_Locate_Point(bounds, p1) < 
                   ST_Line_Locate_Point(bounds, p2);
  
      RAISE DEBUG 'Edge % (left:%, right:%) - ring : % - right_side : %',
        rec.edge_id, rec.left_face, rec.right_face, rrec.path, right_side;

      IF right_side THEN
        right_edges := array_append(right_edges, rec.edge_id);
        old_faceid = rec.right_face;
      ELSE
        left_edges := array_append(left_edges, rec.edge_id);
        old_faceid = rec.left_face;
      END IF;

      IF faceid IS NULL OR faceid = 0 THEN
        faceid = old_faceid;
        old_edgeid = rec.edge_id;
      ELSIF faceid != old_faceid THEN
        RAISE EXCEPTION 'Edge % has face % registered on the side of this face, while edge % has face % on the same side', rec.edge_id, old_faceid, old_edgeid, faceid;
      END IF;

      -- Collect all edges for final full coverage check
      all_edges = ST_Collect(all_edges, rec.geom);

    END LOOP; -- }
  END LOOP; -- }

  IF all_edges IS NULL THEN
    RAISE EXCEPTION 'Found no edges on the polygon boundary';
  END IF;

  RAISE DEBUG 'Left edges: %', left_edges;
  RAISE DEBUG 'Right edges: %', right_edges;

  --
  -- Check that all edges found, taken togheter,
  -- fully match the ring boundary and nothing more
  --
  -- If the test fail either we need to add more edges
  -- from the polygon ring or we need to split
  -- some of the existing ones.
  -- 
  bounds = ST_Boundary(apoly);
  IF NOT ST_isEmpty(ST_SymDifference(bounds, all_edges)) THEN
    IF NOT ST_isEmpty(ST_Difference(bounds, all_edges)) THEN
      RAISE EXCEPTION 'Polygon boundary is not fully defined by existing edges at or near point %', ST_AsText(ST_PointOnSurface(ST_Difference(bounds, all_edges)));
    ELSE
      RAISE EXCEPTION 'Existing edges cover polygon boundary and more at or near point % (invalid topology?)', ST_AsText(ST_PointOnSurface(ST_Difference(all_edges, bounds)));
    END IF;
  END IF;

  IF faceid IS NOT NULL AND faceid != 0 THEN
    IF NOT force_new THEN
      RAISE DEBUG 'Face already known as %, not forcing a new face', faceid;
      RETURN faceid;
    ELSE
      RAISE DEBUG 'Face already known as %, forcing a new face', faceid;
    END IF;
  END IF;

  --
  -- Get new face id from sequence
  --
  FOR rec IN EXECUTE 'SELECT nextval(' ||
    quote_literal(
      quote_ident(atopology) || '.face_face_id_seq'
    ) || ')'
  LOOP
    faceid = rec.nextval;
  END LOOP;

  --
  -- Insert new face 
  --
  EXECUTE 'INSERT INTO '
    || quote_ident(atopology)
    || '.face(face_id, mbr) VALUES('
    -- face_id
    || faceid || ','
    -- minimum bounding rectangle
    || quote_literal(ST_Envelope(apoly)::text)
    || ')';

  --
  -- Update all edges having this face on the left
  --
  IF left_edges IS NOT NULL THEN
    EXECUTE 'UPDATE '
    || quote_ident(atopology)
    || '.edge_data SET left_face = '
    || quote_literal(faceid)
    || ' WHERE edge_id = ANY('
    || quote_literal(left_edges)
    || ') ';
  END IF;

  --
  -- Update all edges having this face on the right
  --
  IF right_edges IS NOT NULL THEN
    EXECUTE 'UPDATE '
    || quote_ident(atopology)
    || '.edge_data SET right_face = '
    || quote_literal(faceid)
    || ' WHERE edge_id = ANY('
    || quote_literal(right_edges)
    || ') ';
  END IF;


  --
  -- Set left_face/right_face of any contained edge 
  --
  EXECUTE 'UPDATE '
    || quote_ident(atopology)
    || '.edge_data SET right_face = '
    || quote_literal(faceid)
    || ', left_face = '
    || quote_literal(faceid)
    || ' WHERE ST_Contains('
    || quote_literal(apoly::text)
    || ', geom)';

  -- 
  -- Set containing_face of any contained node 
  -- 
  EXECUTE 'UPDATE '
    || quote_ident(atopology)
    || '.node SET containing_face = '
    || quote_literal(faceid)
    || ' WHERE containing_face IS NOT NULL AND ST_Contains('
    || quote_literal(apoly::text)
    || ', geom)';

  --
  -- TODO:
  -- Set next_left_face and next_right_face !
  -- These are required by the model, but not really used
  -- by this implementation...
  -- NOTE: should probably be done when adding edges rather than
  --       when registering faces
  --
  RAISE WARNING 'Not updating next_{left,right}_face fields of face boundary edges';


  RETURN faceid;
	
END
$$
LANGUAGE 'plpgsql' VOLATILE;
--} AddFace
