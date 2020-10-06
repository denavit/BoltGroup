%clear all; close all; clc;

x = [0 3 0 3 0 3];
y = [0 0 3 3 6 6];

BG = BoltGroup(x,y);


%%
xP = 1.5+10;
yP = 3;
theta = 0;
[Pn,IC] = BG.Pn_IC(xP,yP,theta)

figure
BG.plot(IC,xP,yP,theta);
axis equal

%%
% d = 1; % counterclockwise positive
% [Mn,IC] = BG.Mn_IC(d)
% 
% figure
% BG.plot(IC,d);
% axis equal
