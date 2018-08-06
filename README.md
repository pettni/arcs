# ARCS

Toolbox in development for abstraction-refinement-based incremental synthesis of correct-by-construction switching protocols.

## Functionality

 - Synthesis methods for augmented finite transition systems.
 - Abstraction tools based on hyper boxes.
 - Demonstrating examples.

## Requirements

 - Matlab (2015b or newer recommended)
 	- Communications Systems Toolbox
 - Yalmip and a supported solver is required for semidefinite optimization, tested version R20160930.
 - [Mosek](https://mosek.com) for conic optimization, tested with version 8.
 - [arrow.m](https://www.mathworks.com/matlabcentral/fileexchange/278-arrow) for certain plotting.
 - Windows 64bit and Linux 64bit for BDD support (so far) or MATLAB mex compiler to compile mex function
 
## Usage

Add the folder ```abstr-ref``` to the Matlab path: ```addpath /path/to/abstr-ref```.

Run tests:
```
cd tests/
runtests
```
Run an example:
```
cd examples_2017/
radiant
```

## Recompile mex files

If for some reason the pre-compiled mex files do not work, run the script ```compile_mex``` in the ```mex_files``` folder. This should create a new mex file that can be moved to the ```abstr-ref/``` folder.

## Authors

Petter Nilsson, Caltech, pettni@caltech.edu

Necmiye Ozay, University of Michigan

Oscar Lindvall Bulancea, KTH

## Todo list
 - Update old examples in `examples/` to work with new software
 - Add custom sos through PolyLinTrans
 - Class for semialgebraic partition/sets (c.f. hyperplane arrangement)
 - Additional tests
	- Dual algorithms
    - Control strategies
    - Candidate sets
    - System-level test
 - Improve documentation
 - Write custom SOS optimization to avoid overhead (translate python socp impl.)
 - Improve documentation and type security of mex code

## References

 - P Nilsson and N Ozay, "Incremental Synthesis of Switching Protocols via Abstraction Refinement", Proc. of IEEE CDC, 2014.
 - P Nilsson, N Ozay and J Liu, "Augmented Finite Transition Systems as Abstractions for Control Synthesis", Journal of Discrete Event Dynamic Systems (Special issue on Formal Methods in Control), 2017.
 - O B Lindvall, P Nilsson and N Ozay, "Nonuniform abstractions, refinement and controller synthesis with novel BDD encodings", in Proceedings of the IFAC Conference on Analysis and Design of Hybrid Systems, 2018 
 
 ## Acknowledgments: 
 This research is supported in part by DARPA grant N66001-14-1-4045.
