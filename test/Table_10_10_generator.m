clear all; close all; clc; 

list_Fy = [36 50];
list_d  = [3/4 7/8 1 1+1/8];
list_n  = 12:-1:2;
list_bg = {'A','B'};
list_tc = {'N','X'};
list_tp = [1/4 5/16 3/8 7/16 1/2 9/16];

list_edge_dist = [1 1+1/8 1+1/4 1+1/2];

fid = fopen('Table_10-10.csv','w');
fprintf(fid,'Fy,d,n,bolt_group,thread_cond,tp,fRn,limit_state\n');

for iFy = 1:length(list_Fy)
    for id = 1:length(list_d)
        for in = 1:length(list_n)
            Fy = list_Fy(iFy);
            d  = list_d(id);
            n  = list_n(in);
            fprintf('Working on Fy = %g ksi, d = %g in, n = %i\n',Fy,d,n);
            for ibg = 1:length(list_bg)
                for itc = 1:length(list_tc)
                    for itp = 1:length(list_tp)
                        bolt_group  = list_bg{ibg};
                        thread_cond = list_tc{itc};
                        hole_type   = 'STD';
                        tp = list_tp(itp);

                        conn = SinglePlateShearConnection(Fy,d,n,bolt_group,thread_cond,hole_type,tp);
                        conn.use_tabulated_C = true;
                        
                        fprintf(fid,'%g,%g,%i,%s,%s,%g,',Fy,d,n,bolt_group,thread_cond,tp);
                        if tp > conn.max_tp
                            fprintf(fid,'---,---\n');
                        else
                            assert(conn.is_conventional(),'Connection must be conventional');
                            
                            % conn.bolt_strength_method_override = 'Bearing and Tearout as Concentric';
                            [fRn,limit_state] = conn.R('design');

                            fprintf(fid,'%g,%s\n',fRn,limit_state);
                        end
                        
                    end
                end
            end
        end
    end
end

fclose(fid);
