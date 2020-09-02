# ORC
Optimised Relaxed Clock


## Installation instructions

### BEAUTI



### By hand

This assumes you have already installed BEAST2 on your machine
https://www.beast2.org/

To install this package manually, first clone this repository and all of its dependencies from GitHub, including BEAST2

```
mkdir beast2
cd beast2
git clone https://github.com/CompEvol/BEAST2
git clone https://github.com/BEAST2-Dev/BEASTLabs
git clone https://github.com/Rong419/ConstantDistanceOperator
git clone https://github.com/jordandouglas/ORC
```

Then navigate to into each directory and install them using ant

```
cd BEAST2
ant linux

cd BEASTLabs
ant addon

cd ConstantDistanceOperator
ant addon

cd ORC
ant addon

cd ..
```

Finally, extract the contents of these addons:

```
unzip -o BEASTLabs/build/dist/ORC.addon.v*zip -d ~/.beast/2.6/BEASTLabs/.
unzip -o ConstantDistanceOperator/build/dist/ORC.addon.v*zip -d ~/.beast/2.6/FastRelaxedClockLogNormal/.
unzip -o ORC/build/dist/ORC.addon.v*zip -d ~/.beast/2.6/ORC/.
```


This assumes that BEAST2 is installed in the ~/.beast2/ directory.

If you are not using BEAST 2.6 then replace 2.6 with the appropriate version number.
For more details see https://beast2.blogs.auckland.ac.nz/managing-packages/#Install_by_hand .


## Automated generation of narrow exchange operators
Run matlab/GenerateOperators.m using MatLab to generate all 53 / 63 solvable operators 





