clear all; close all; clc;

%% Design Example II.A-17A
% Design Strength 
%  52.2 k  [Table 10.10a]

data.Fy          = 36;
data.Fu          = 58;
data.d           = 0.75;
data.n           = 4;
data.bolt_group  = 'A';
data.thread_cond = 'N';
data.hole_type   = 'STD';
data.tp          = 0.25;
data.a           = 3;
data.s           = 3;
data.g           = 3;
data.nr          = 1;
data.leh         = 1.5;
data.lev         = 1.25;

conn = SinglePlateShearConnection(data);
% conn.check_all_limit_states_for_conventional = true;

[fRn1,limit_state1,results1] = conn.R('design');

conn.plot_IC_results();

for j = 1:length(results1.strengths)
    fprintf('%8.02f kips (%s)\n',results1.strengths(j),results1.limit_states{j});
end