![alt text](https://github.com/rdzhao/ESES/blob/master/fig/pipeline.png)

# ESES
Software for Eulerian [solvent excluded surface](https://www.annualreviews.org/doi/abs/10.1146/annurev.bb.06.060177.001055) (SES), which geneates the SES in a Eulerian regular grid with parallel treatment.

# Compilation
Only OpenMP is required for parallel computation. 
Make sure OpenMP is installed before compiling the source code.

#Usage
After compilation, An application "MS_Intersection" is generated. 
User should provide 4 parameters.
1) [.pqr file](https://www.mdanalysis.org/docs/documentation_pages/coordinates/PQR.html) decrsbing the molecules
2) probe size
3) resolution
4) grid margin

#Output
There are 3 output files.
1) bounding_box.txt: describing basic grid layout information.
2) intersection_info.txt: describing intersctions by providing which grid edge the intersection lie on, the intersection point, normal, and which type of surface patch define by SES it belongs to.
3) grid_info.txt: descrbing whether the grid points are inside(1) or outside(-1);


