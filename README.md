# BoltGroup
This project contains MATLAB functions for determining the strength of eccentrically loaded bolt groups using the instantaneous center of rotation method. This code was initially developed as part of a CE571 term project by Walker Trent.

The main calculations are implemented in the BoltGroup.m MATLAB class. Functions in this class use [fsolve](https://www.mathworks.com/help/optim/ug/fsolve.html) (part of the Optimization Toolbox) to iteratively determine the location of the instantaneous center of rotation.
