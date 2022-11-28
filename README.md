About
-----

The codes here are run in dokcer (image: underworldcode/underworld2:2.13.0b)

- RefM1_noUniformmesh: 
The problems of the shape.Polygon in non uniform mesh that some swarms are not assigned to the material index.

- 1_23_02_FreeSurface_Kaus2010_Rayleigh-Taylor_Instability_NUMesh: 
Swam number increases a lot in non-uniform mesh.
The value for each particlesPerCell is 40. In the uniform mesh, the swarm numbers are 100k (250*40) initially and increase to 105k when the model ends running,  the swarm number in cell will be 40+-5 when the model ends;
But for the non-uniform mesh (the half part of the model mesh is sparse), the swarm numbers are 100k initially and increase to 130k when the model ends running,  the swarm number in some cells will be 100+



