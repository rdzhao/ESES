# ESES
Software for Eulerian [solvent excluded surface](https://www.annualreviews.org/doi/abs/10.1146/annurev.bb.06.060177.001055) (SES), which geneates the SES in a Eulerian regular grid with parallel treatment.

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

# Output
There are 3 output files.
1) bounding_box.txt: describing basic grid layout information.
2) intersection_info.txt: describing intersctions by providing which grid edge the intersection lie on, the intersection point, normal, and which type of surface patch define by SES it belongs to.
3) grid_info.txt: descrbing whether the grid points are inside(1) or outside(-1);
