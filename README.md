# Radiative transfer code Dusty

[Team space:](https://github.com/dirac-institute/DustyCollaboration)


## There are two versions: 

1) The original 1999 version, known internally as V2. Instructions,
   the code, and supporting files for V2 are available from
   [the old Dusty's website](http://faculty.washington.edu/ivezic/dusty_web/)

2) The latest and greatest ("new and improved") version, known internally as V4. 
   You can get V4 as *tar  
   [file distribution](release/dusty.tar), 
   which also includes [the V4 manual](release/dusty/docs/manual.pdf)
   
## Important note:

V4 is not nearly as much tested and verified as V2. V2 is recommended over V4, unless
you require V4 functionality that doesn't exist in V2:

- slab geometry (V2 only supports spherical geometry)
- external heating (V2 only supports a central source of radiation)
- multi-grain mixtures (V2 uses an approximation that is still single grain type model) 


## Building steps for V4: 

i) To generate V4 release files (directory release), run 
> sh generate_release.sh

ii) then compile
> cd release/dusty

> make

NB: for MacOS users - if you encounter "ld: library not found for -lcrt1.o" problem, 
then do
> xcode-select --install

iii) and finally run the example master file:
> ./dusty dusty.mas


