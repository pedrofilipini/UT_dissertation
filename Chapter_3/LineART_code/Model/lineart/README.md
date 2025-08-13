# multibart R Package

## About
Implements a variety of BART models using a regression on basis functions in the leaves.

Function bcf_new is for TS-BART.
Function lbart_new is for Linear-BART.

This is a version of lbart_new where each leaf has a simple linear regression by using the same framework of basis functions.

This version of lbart_new works (functions for grow and prune have been modified, as well as metropolis, and fit). The m_is are counting if the variables have been used on each tree.

Still have to deal with the initialization of the parameters of the leaves (MLE) on lbart_new.

The prediction function for the lbart model works, but I still need to include locks to avoid the user from "breaking" the function (like adding a matrix of the wrong size inside the function).

Since it has been adapted from the BART package, the prediction of the lbart should still work in parallel when used on Linux (I have not tested this yet, but it should work).

## Documentation

## Installation

## Author
