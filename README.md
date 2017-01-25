# ABSTR-REFINEMENT

Toolbox in development for incremental synthesis of correct-by-construction switching protocols. 

## Functionality

 - Synthesis methods for augmented finite transition systems.
 - Abstraction tools based on hyper boxes.
 - Demonstrating examples.

## Requirements

 - Matlab 2015b or newer. Older versions work if the calls to ```builtin('_ismemberhelper',aPost,X)``` in ```pre.m``` are replaced with ```ismember(aPost, X)```, but this will be slower.
 - Yalmip is required for semidefinite optimization, tested version R20160930.

## Usage

Add the folder ```abstr-ref``` to the Matlab path: ```addpath /path/to/abstr-ref```.

Run tests:
```
run(tests)
```
Run an example:
```
cd examples/
radiant_new
```

## Authors

Petter Nilsson, University of Michigan, pettni@umich.edu

Necmiye Ozay, University of Michigan

## Todo list

 - Implement remaining algorithms from manuscript, including control extraction
 - Create class that handles control strategy
 - (Candidate) losing set algorithms
 - Recursive candidate set computation
 - Class for semialgebraic partition/sets
 - Additional tests
    - Control strategies
    - System-level test

## References

 - Nilsson, P and Ozay, N, "Incremental Synthesis of Switching Protocols via Abstraction Refinement", Proc. of IEEE CDC, 2014,
