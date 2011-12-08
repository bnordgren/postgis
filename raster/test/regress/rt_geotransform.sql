SELECT (gt).imag as imag, 
       (gt).jmag as jmag, 
       degrees((gt).theta_i) as theta_i, 
       degrees((gt).theta_ij) as theta_ij,
       (gt).xoffset as xoffset, 
       (gt).yoffset as yoffset 
   FROM (SELECT ST_GetGeotransform(ST_MakeEmptyRaster(10,10, 0,0, 20)) as gt) as dummy  ; 

SELECT (gt).imag as imag, 
       (gt).jmag as jmag, 
       degrees((gt).theta_i) as theta_i, 
       degrees((gt).theta_ij) as theta_ij,
       (gt).xoffset as xoffset, 
       (gt).yoffset as yoffset 
   FROM 
      (SELECT 
          ST_GetGeotransform(
             ST_SetScale(
               ST_SetSkew(
                 ST_MakeEmptyRaster(10,10,0,0,20), 
                 -sqrt(2), sqrt(2)), 
               sqrt(2))) as gt) as dummy ;


SELECT (gt).imag as imag, 
       (gt).jmag as jmag, 
       degrees((gt).theta_i) as theta_i, 
       degrees((gt).theta_ij) as theta_ij,
       (gt).xoffset as xoffset, 
       (gt).yoffset as yoffset 
  FROM (SELECT ST_GetGeotransform(
               ST_SetGeotransform(ST_MakeEmptyRaster(10,10, 0,0, 20), 
                                  50, 85, radians(45), radians(90), 20, 30)) as gt) as foo ; 

                                  
                                  