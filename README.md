# ESES
Software for Eulerian [solvent excluded surface](https://www.annualreviews.org/doi/abs/10.1146/annurev.bb.06.060177.001055) (SES), which geneates the SES in a Eulerian regular grid with parallel treatment.
Our software support Catesian/Lagrangian mesh generation with files for visualization.

![alt text](https://github.com/rdzhao/ESES/blob/master/fig/pipeline.png)
![alt text](https://github.com/rdzhao/ESES/blob/master/fig/huge.png)

# Compilation
The code is tested only on MacOS/Linux, but it should require only minimal work to get it working on Windows.
Only OpenMP is required as an external dependency. 
Before compiling, you need to install OpenMP.
#### Install OpenMP in MacOS 
```
brew install llvm
brew install libomp
```
#### Install OpenMP in Linux
```
sudo apt-get install libomp-dev
```
#### Compilation procedure
```
cd ESES_root
mkdir build
cd ./build
cmake ..
make
```

# Usage
After compilation, An application "MS_Intersection" is generated. 
User should provide 4 parameters.
1) [.pqr file](https://www.mdanalysis.org/docs/documentation_pages/coordinates/PQR.html) decrsbing the molecules
2) probe size (in angstrom)
3) resolution (in angstrom)
4) grid margin (in angstrom)

For example, if we are given a .prq file "1ajj.pqr", the command looks like
```
./MS_Intersection 1ajj.pqr 1.4 0.4 1.0
```

# Output files
The output file includes Cartesian/Langrangian meshes with files for visualization.
#### Cartesian Mesh
**bounding_box.txt**: describing basic grid layout information.
```
x_min y_min z_min       # min corner
x_max y_max z_max       # max corner
x_dim y_dim z_dim       # grid dimension
```
**intersection_info.txt**: describing intersctions by providing which grid edge the intersection lie on, the intersection point, normal, and which type of surface patch define by SES it belongs to.
**gp1** and **pg2** are two grid points describing the grid edge.
**p** is the intersection position vector.
**n** is the normal vector.
**atom_k** are adjacent atoms of the intersection. Adjacent atoms are a descripiion of which type of patches (convex, saddel, concave) the intersection lies on. Note that convex points have 1 adjacent atom, saddle points have 2 adjacent atoms, and concave points have 3 adjacent atoms. 

```
gp1_x gp1_y gp1_z gp1_x gp1_y gp1_z p_x p_y p_z n_x n_y n_z atom_0 atom_1 
...
```
**grid_info.txt**: descrbing whether the grid points are inside(1) or outside(-1).
**pg** is the grip point.
```
gp_x gp_y gp_z inside/outside
...
```

#### Lagrangian Mesh
**marching_cubes.obj**: the Lagrangian mesh is saved as a .obj file.

#### files for visualziation
We also have all the related files for visualizing space fill of atoms, different type of patches, cartesian mesh in grid and partition of surface by different chain. 
