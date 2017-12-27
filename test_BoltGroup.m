clear all; close all; clc;

x = [0 3 0 3 0 3];
y = [0 0 3 3 6 6];

BG = BoltGroup(x,y);


%%
[Pn,IC] = BG.Pn_IC(4,3,0)

figure
BG.plot(IC);
axis equal

%%
[Mn,IC] = BG.Mn_IC()

figure
BG.plot(IC);
axis equal
