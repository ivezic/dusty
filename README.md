# Radiative transfer code Dusty

There are two versions: 

1) The original 1999 version, known internally as V2. Instructions,
   the code, and supporting files for V2 are available from
   [the old Dusty's website](http://faculty.washington.edu/ivezic/dusty_web/)

2) The latest and greatest ("new and improved") version, known internally as V4. 
   You can get V4 as *tar  
   [file distribution](release/dusty.tar), 
   which also includes the manual as [docs/manual.pdf](docs/manual.pdf)


Internal info: 
i) To generate V4 release files in directory release/, run 
> sh generate_release.sh

ii) then compile
> cd release/dusty
> make

(for MacOS users: if you encounter "ld: library not found for -lcrt1.o" problem, 
then do
> xcode-select --install) 

iii) and finally run the example master file:
./dusty dusty.mas


