-----------------------------------------------------------------------
-- $Id$
--
-- Copyright (c) 2009 Sandro Santilli <strk@keybit.net>
--
-- This is free software; you can redistribute and/or modify it under
-- the terms of the GNU General Public Licence. See the COPYING file.
-----------------------------------------------------------------------

CREATE TABLE rt_box2d_test (
    id numeric,
    name text,
    rast raster,
    env box2d
);

-- 10x20, ip:0.5,0.5 scale:2,3
INSERT INTO rt_box2d_test 
VALUES ( 0, '10x20, ip:0.5,0.5 scale:2,3 skew:0,0',
(
'01' -- little endian (uint8 ndr)
|| 
'0000' -- version (uint16 0)
||
'0000' -- nBands (uint16 0)
||
'0000000000000040' -- scaleX (float64 2)
||
'0000000000000840' -- scaleY (float64 3)
||
'000000000000E03F' -- ipX (float64 0.5)
||
'000000000000E03F' -- ipY (float64 0.5)
||
'0000000000000000' -- skewX (float64 0)
||
'0000000000000000' -- skewY (float64 0)
||
'0A000000' -- SRID (int32 10)
||
'0A00' -- width (uint16 10)
||
'1400' -- height (uint16 20)
)::raster
,'BOX(0.5 0.5,20.5 60.5)' -- expected envelope (20x60) == (10*2 x 20*3)
);

INSERT INTO rt_box2d_test 
VALUES ( 1, '1x1, ip:2.5,2.5 scale:5,5 skew:0,0',
(
'01' -- little endian (uint8 ndr)
|| 
'0000' -- version (uint16 0)
||
'0000' -- nBands (uint16 0)
||
'0000000000001440' -- scaleX (float64 5)
||
'0000000000001440' -- scaleY (float64 5)
||
'0000000000000440' -- ipX (float64 2.5)
||
'0000000000000440' -- ipY (float64 2.5)
||
'0000000000000000' -- skewX (float64 0)
||
'0000000000000000' -- skewY (float64 0)
||
'00000000' -- SRID (int32 0)
||
'0100' -- width (uint16 1)
||
'0100' -- height (uint16 1)
)::raster
,'BOX(2.5 2.5,7.5 7.5)' -- expected envelope 
);

INSERT INTO rt_box2d_test 
VALUES ( 2, '1x1, ip:7.5,2.5 scale:5,5 skew:0,0',
(
'01' -- little endian (uint8 ndr)
|| 
'0000' -- version (uint16 0)
||
'0000' -- nBands (uint16 0)
||
'0000000000001440' -- scaleX (float64 5)
||
'0000000000001440' -- scaleY (float64 5)
||
'0000000000001E40' -- ipX (float64 7.5)
||
'0000000000000440' -- ipY (float64 2.5)
||
'0000000000000000' -- skewX (float64 0)
||
'0000000000000000' -- skewY (float64 0)
||
'00000000' -- SRID (int32 0)
||
'0100' -- width (uint16 1)
||
'0100' -- height (uint16 1)
)::raster
,'BOX(7.5 2.5,12.5 7.5)' -- expected envelope 
);
