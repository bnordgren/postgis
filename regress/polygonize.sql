SELECT 1, ST_AsText(ST_Polygonize(
'GEOMETRYCOLLECTION(
  MULTILINESTRING(
    (1656318.45 4833344.45,1656321.79 4833339.62,1656312.54 4833333.49),
    (1656312.54 4833333.49,1656309.68 4833337.07)
  ),
  LINESTRING(1656309.68 4833337.07,1656318.45 4833344.45)
)'::geometry));
