## License

Copyright (c) 2021, Tilo Arens, Roland Griesmaier, Marvin Knöller

This software is released under GNU General Public License, Version 3.
The full license text is provided in the file LICENSE included in this repository 
or can be obtained from http://www.gnu.org/licenses/


## Content of this project
This is a guide to generate the figures that have been used in the work

**„Maximizing the electromagnetic chirality of thin dielectric tubes“**

_By Tilo Ares, Roland Griesmaier and Marvin Knöller._

You find all needed Matlab files in the folders. 

## An overview:
- [] The folder **AuxiliaryUtilities** includes the files _dotReal_ (dot-product)
								   _eval_phi_ (needed for the BFGS method)
								   _Tensor_Mat_Mult_ (multiplies 2 3-tensors slice by slice along the third dimension)
- [] The folder **ChiralityUtilities** includes the files _chiral_ (computes chirality measures)
								  _FarFieldMatrixFunction_(s) (computes FarFieldOperator approximations given parameters and a curve)
- [] The folder **Derivatives** includes all derivates that are needed in the BFGS method especially derivates of 
-far field operators
-all penalty terms
-entire wave fields
-chirality measures
- []The folder **PaperPlots** includes plotting routines
- []The folder **Penalties** includes the penalty terms given a spline
- []The folder **SplineUtilities** includes all routines to generate cubic splines in the way they are needed for the other programs 
——

——
For generating the figures 6.1 - 6.6 from the work, run

generateFigure6_1
generateFigure6_2
generateFigure6_3
generateFigure6_4
generateFigure6_5 and
generateFigure6_6.

If you do not wish to run the optimization from the start, you can set the value generateFiles to 0.
This will generate the figures from pre-computed files.
——

——
The computations have been carried out on a MacBook Pro 15“ model from 2018.
Generating all figures from scratch took approximately 2-3 days.
Computations have been carried out using the Matlab 2020a version.

The code uses parallelization from the Matlab Parallelization Toolbox.
——
