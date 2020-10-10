clear all; close all; clc;

%% Design Example II.A-19A
% Design Strength 
%  41.7 k  [Poison Bolt Method, Bolt Shear Rupture]
% 130 k    [Shear Yielding of Plate]
% 111 k    [Shear Rupture of Plate]
% 117 k    [Block Shear Rupture of Plate]
%  57.9 k  [Plate Interaction, (36 k)/sqrt(0.386)]
%  64.8 k  [Lateral-Torsional Buckling, (583 k-in.)/(9 in.)]
%  61.9 k  [Flexural Rupture of Plate, (557 k-in.)/(9 in.)]

data.Fy          = 36;
data.Fu          = 58;
data.d           = 0.75;
data.n           = 4;
data.bolt_group  = 'A';
data.thread_cond = 'N';
data.hole_type   = 'STD';
data.tp          = 0.5;
data.a           = 9;
data.s           = 3;
data.g           = 3;
data.nr          = 2;
data.leh         = 1.25;
data.lev         = 1.5;

conn = SinglePlateShearConnection(data);

conn.bolt_strength_method_override = 'Poison Bolt';
[fRn1,limit_state1,results1] = conn.R('design');

conn.plot_IC_results();

for j = 1:length(results1.strengths)
    fprintf('%8.02f kips (%s)\n',results1.strengths(j),results1.limit_states{j});
end