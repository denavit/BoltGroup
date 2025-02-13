clear all; close all; clc;

% Bolt Group 1 (compare to results in Figure 4) 
%
% Reynolds, M., Redl, E., and Uang, C.-M. (2021). “Effect of Bolt-Hole 
% Clearance on the Ultimate Strength of Eccentrically Loaded Bolt Groups.” 
% Journal of Structural Engineering, ASCE, 147(1), 04020285.

x               = [0 0 0 0 0];
y               = [-2 -1 0 1 2]*3;
Rult            = 1;

BG = BoltGroup(x,y,Rult);
[Pn0,IC0] = BG.Pn_IC(-12,0,0);
BG.plot(IC0,-12,0,0)

BG.load_deformation_type = 'standard_with_slip';
BG.Rslip = 0;
BG.Dslip = 1/8;
[PnS,ICS] = BG.Pn_IC(-12,0,0);
BG.plot(ICS,-12,0,0)