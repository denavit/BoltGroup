classdef SinglePlateShearConnection
    
    properties
        Fy          % Plate yield strength
        Fu          % Plate ultimate strength
        d           % Bolt diameter
        n           % Number of bolts in one vertical row
        bolt_group
        thread_cond
        hole_type   % Currently only impemented for STD holes
        tp          % Plate thickness
        a           % Distance from bolt line to weld line
        s           % Spacing of bolts
        g           % Horizontal spacing between vertical rows of bolts
        nr          % Number of vertical rows of bolts
        leh
        lev
        
        
        % Calculation Parameters
        default_strength_type                   = 'nominal';
        conventional_override                   = [];
        check_all_limit_states_for_conventional = false;
        bolt_strength_method_override           = [];
        use_tabulated_C                         = false;
        eccentricity_on_right_side              = true;
        bolt_hole_def_considered                = false;
    end
    
    methods
        function obj = SinglePlateShearConnection(Fy,d,n,bolt_group,thread_cond,hole_type,tp)
            if nargin == 1
                % Create a SinglePlateShearConnection with general
                % parameters defined with a structure
                data = Fy;
                obj.Fy          = data.Fy;
                obj.Fu          = data.Fu;
                obj.d           = data.d;
                obj.n           = data.n;
                obj.bolt_group  = data.bolt_group;
                obj.thread_cond = data.thread_cond;
                obj.hole_type   = data.hole_type;
                obj.tp          = data.tp;
                obj.a           = data.a;
                obj.s           = data.s;
                obj.g           = data.g;
                obj.nr          = data.nr;
                obj.leh         = data.leh;
                obj.lev         = data.lev;
            else
                % Create a SinglePlateShearConnection based on the Table
                % 10-10 configuration 
                obj.Fy          = Fy;
                
                % Set Fu based on Fy
                if Fy == 36
                    obj.Fu = 58;
                elseif Fy == 50
                    obj.Fu = 65;
                else
                    error('Unknown Fy: %g',Fy)
                end                
                
                obj.d           = d;
                obj.n           = n;
                obj.bolt_group  = bolt_group;
                obj.thread_cond = thread_cond;
                obj.hole_type   = hole_type;
                obj.tp          = tp;

                obj.a           = 3;
                obj.s           = 3;
                obj.g           = 3;
                obj.nr          = 1;
                
                % Set leh based on d
                obj.leh = 2*d;

                % Set lev based on d
                if d == 0.75
                    obj.lev = 1.25;
                elseif d == 0.875
                    obj.lev = 1.5;
                elseif d == 1
                    obj.lev = 1.75;
                elseif d == 1.125
                    obj.lev = 2.0;
                else
                    error('Unknown d: %g',d)
                end
            end
        end
        function d = dp(obj) % depth of plate
            d = (obj.n-1)*obj.s + 2*obj.lev;         
        end
        function tf = is_conventional(obj)
            if ~isempty(obj.conventional_override)
                tf = obj.conventional_override;
            else
                tf1 = obj.nr == 1;
                tf2 = (obj.n >= 2) && (obj.n <= 12);
                tf3 = obj.a <= 3.5;
                tf4 = strcmpi(obj.hole_type,'STD') || strcmpi(obj.hole_type,'SSLT');
                tf5 = obj.lev >= minimum_edge_distance(obj.d,obj.hole_type);
                tf6 = obj.leh >= 2*obj.d;
                tf7 = obj.tp <= obj.max_tp;
                tf  = all([tf1 tf2 tf3 tf4 tf5 tf6 tf7]);
            end
        end
        function t = max_tp(obj)
            switch obj.hole_type
                case 'STD'
                    if obj.n <= 5
                        t = obj.d/2 + 1/16;
                    else
                        t = obj.d/2 - 1/16;
                    end
                case 'SSLT'
                    if obj.n <= 5
                        t = Inf;
                    else
                        t = obj.d/2 + 1/16;
                    end
                otherwise
                    error('Unknown hole type: %s',obj.hole_type);
            end
        end
        function ecc = e(obj)
            % Eccentricity from first row of bolts
            if obj.is_conventional
                switch obj.hole_type
                    case 'STD'
                        if obj.n <= 5
                            ecc = obj.a/2;
                        else
                            ecc = obj.a;
                        end 
                    case 'SSLT'
                        ecc = obj.a/2;
                    otherwise
                        error('Unknown hole type: %s',obj.hole_type);
                end
            else
                ecc = obj.a;
            end
            
            if obj.eccentricity_on_right_side
                ecc = ecc + obj.g*(obj.nr-1);
            else
                ecc = -ecc;
            end         
        end
        function a = Ab(obj)
            a = pi/4*obj.d^2;
        end
        function dh = dh(obj)
            if obj.d <= 7/8
                dh = obj.d + 1/16;
            else
                dh = obj.d + 1/8;
            end
        end
        function z = Znet(obj)
            z = obj.tp*obj.dp^2/4;
            y = ((0:(obj.n-1))-((obj.n-1)/2))*obj.s;
            for i = 1:obj.n
                if y(i) == 0
                    z = z - obj.tp*(obj.dh+1/16)^2/4;
                else
                    z = z - obj.tp*(obj.dh+1/16)*abs(y(i));
                end
            end
        end
            
        function [x,y] = get_bolt_positions(obj)
            x = nan(1,obj.n*obj.nr);
            y = nan(1,obj.n*obj.nr);
            for i = 1:obj.nr
                ind = (1:obj.n) + (i-1)*obj.n;
                x(ind) = obj.g*(i-1)*ones(1,obj.n);
                y(ind) = ((0:(obj.n-1))-((obj.n-1)/2))*obj.s;
            end
        end
        function fnv = Fnv(obj)
            switch obj.bolt_group
                case 'A'
                    switch obj.thread_cond
                        case 'N'; fnv = 54;
                        case 'X'; fnv = 68;
                        otherwise; error('Unknown thread condition: %s',obj.thread_cond);
                    end
                case 'B'
                    switch obj.thread_cond
                        case 'N'; fnv = 68;
                        case 'X'; fnv = 84;
                        otherwise; error('Unknown thread condition: %s',obj.thread_cond);
                    end
                otherwise
                    error('Unknown thread conddition: %s',obj.thread_cond);
            end
        end
        function [controlling_R,controlling_limit_state,results] = R(obj,strength_type)
            if nargin < 2
                strength_type = obj.default_strength_type;
            end
            strengths    = zeros(0,1);
            limit_states = cell(0,1);
            
            if isempty(obj.bolt_strength_method_override)
                % Defaults
                if obj.is_conventional
                    bolt_strength_method = 'Bearing and Tearout as Concentric';
                else
                    bolt_strength_method = 'Modified IC';
                end
            else
                bolt_strength_method = obj.bolt_strength_method_override;
            end
            
            switch bolt_strength_method
                case 'Bearing and Tearout as Concentric'

                    % Bolt shear (as eccentric) 
                    if obj.use_tabulated_C
                        C   = C_from_table(obj.n);
                        Rn  = obj.Fnv*obj.Ab*C;
                    else
                        [x,y] = obj.get_bolt_positions;
                        Rult = obj.Fnv*obj.Ab;
                        BG  = BoltGroup(x,y,Rult);
                        Rn  = BG.Pn_IC(obj.e,0,0);
                        % fRn/Rult; % this is C
                    end

                    switch lower(strength_type)
                        case 'nominal';     R = Rn;
                        case 'design';      R = 0.75*Rn;
                        case 'available';   R = Rn/2.00;
                        otherwise; error('Unknown strength type: %s',strength_type)
                    end

                    strengths    = vertcat(strengths,R);
                    limit_states = vertcat(limit_states,'Bolt shear (as eccentric)');   

                    % Bearing and Tearout (as concentric)
                    lc  = obj.s - obj.dh;
                    if obj.bolt_hole_def_considered
                        rn1 = min([2.4*obj.d*obj.tp*obj.Fu 1.2*lc*obj.tp*obj.Fu]);
                    else
                        rn1 = min([3.0*obj.d*obj.tp*obj.Fu 1.5*lc*obj.tp*obj.Fu]);
                    end
                    lc  = obj.lev - obj.dh/2;
                    if obj.bolt_hole_def_considered
                        rn2 = min([2.4*obj.d*obj.tp*obj.Fu 1.2*lc*obj.tp*obj.Fu]);
                    else
                        rn2 = min([3.0*obj.d*obj.tp*obj.Fu 1.5*lc*obj.tp*obj.Fu]);
                    end
                    Rn  = (rn1*(obj.n-1) + rn2);

                    switch lower(strength_type)
                        case 'nominal';     R = Rn;
                        case 'design';      R = 0.75*Rn;
                        case 'available';   R = Rn/2.00;
                        otherwise; error('Unknown strength type: %s',strength_type)
                    end
                    
                    strengths    = vertcat(strengths,R);
                    limit_states = vertcat(limit_states,'Bearing and Tearout (as concentric)');   
                
                case 'Modified IC'
                    % Bolt shear, bearing, and tearout (as eccentric) 
                    [x,y] = obj.get_bolt_positions;
                    if obj.bolt_hole_def_considered
                        Rult = min([obj.Fnv*obj.Ab 2.4*obj.d*obj.tp*obj.Fu]);
                        Rult_over_lc = 1.2*obj.tp*obj.Fu;
                    else
                        Rult = min([obj.Fnv*obj.Ab 3.0*obj.d*obj.tp*obj.Fu]);
                        Rult_over_lc = 1.5*obj.tp*obj.Fu;
                    end
                    BG = BoltGroupWithTearout(x,y,Rult,Rult_over_lc,obj.dh);
                    BG.edge_right_x  = (obj.nr-1)*obj.g + obj.leh;
                    BG.edge_top_y    = obj.dp/2;
                    BG.edge_bottom_y = -obj.dp/2;
                    Rn = BG.Pn_IC(obj.e,0,0);
                    % fRn/Rult; % this is C

                    switch lower(strength_type)
                        case 'nominal';     R = Rn;
                        case 'design';      R = 0.75*Rn;
                        case 'available';   R = Rn/2.00;
                        otherwise; error('Unknown strength type: %s',strength_type)
                    end
                    
                    strengths    = vertcat(strengths,R);
                    limit_states = vertcat(limit_states,'Modified IC'); 
                    
                case 'Poison Bolt'
                    lc_min = min([obj.lev-0.5*obj.dh obj.leh-0.5*obj.dh obj.s-obj.dh obj.g-obj.dh]);
                    if obj.bolt_hole_def_considered
                        Rult = min([obj.Fnv*obj.Ab 2.4*obj.d*obj.tp*obj.Fu 1.2*lc_min*obj.tp*obj.Fu]);
                    else
                        Rult = min([obj.Fnv*obj.Ab 3.0*obj.d*obj.tp*obj.Fu 1.5*lc_min*obj.tp*obj.Fu]);
                    end
                    
                    if obj.use_tabulated_C
                        Rn = Rult*C_from_table(obj.n);
                    else
                        [x,y] = obj.get_bolt_positions;
                        BG = BoltGroup(x,y,Rult);
                        Rn = BG.Pn_IC(obj.e,0,0);
                        %Rn/Rult % this is C
                    end

                    switch lower(strength_type)
                        case 'nominal';     R = Rn;
                        case 'design';      R = 0.75*Rn;
                        case 'available';   R = Rn/2.00;
                        otherwise; error('Unknown strength type: %s',strength_type)
                    end

                    strengths    = vertcat(strengths,R);
                    limit_states = vertcat(limit_states,'Poison Bolt Method'); 

                otherwise
                    error('Unknown bolt_strength_method: %s',bolt_strength_method)
            end
            
            % Shear Yielding of the Plate
            Agv = obj.dp*obj.tp;
            Rn  = 0.6*obj.Fy*Agv;
            
            switch lower(strength_type)
                case 'nominal';     R = Rn;
                case 'design';      R = 1.0*Rn;
                case 'available';   R = Rn/1.50;
                otherwise; error('Unknown strength type: %s',strength_type)
            end

            strengths    = vertcat(strengths,R);
            limit_states = vertcat(limit_states,'Shear Yielding of the Plate');             
            
            % Shear Rupture of the Plate
            Anv = (obj.dp-obj.n*(obj.dh+(1/16)))*obj.tp;
            Rn  = 0.6*obj.Fu*Anv;
            
            switch lower(strength_type)
                case 'nominal';     R = Rn;
                case 'design';      R = 0.75*Rn;
                case 'available';   R = Rn/2.00;
                otherwise; error('Unknown strength type: %s',strength_type)
            end

            strengths    = vertcat(strengths,R);
            limit_states = vertcat(limit_states,'Shear Rupture of the Plate');             
            
            % Block Shear of the Plate
            Agv = ((obj.n-1)*obj.s + obj.lev)*obj.tp;
            Anv = Agv - (obj.n-0.5)*(obj.dh+(1/16))*obj.tp;
            Ant = ((obj.nr-1)*obj.g + obj.leh - (obj.nr-0.5)*(obj.dh+(1/16)))*obj.tp;
            if obj.nr == 1
                Ubs = 1.0;
            else
                Ubs = 0.5;
            end
            Rn  = min([0.6*obj.Fu*Anv+Ubs*obj.Fu*Ant 0.6*obj.Fy*Agv+Ubs*obj.Fu*Ant]);
            
            switch lower(strength_type)
                case 'nominal';     R = Rn;
                case 'design';      R = 0.75*Rn;
                case 'available';   R = Rn/2.00;
                otherwise; error('Unknown strength type: %s',strength_type)
            end

            strengths    = vertcat(strengths,R);
            limit_states = vertcat(limit_states,'Block Shear Rupture of the Plate'); 
            
            if ~obj.is_conventional || obj.check_all_limit_states_for_conventional
                % Interaction Strength of the Plate
                Zg  = obj.tp*obj.dp^2/4;
                Mn  = obj.Fy*Zg;
                Agv = obj.dp*obj.tp;
                Vn  = 0.6*obj.Fy*Agv;

                switch lower(strength_type)
                    case 'nominal';     R = 1/sqrt((1/Vn)^2+(obj.a/Mn)^2);
                    case 'design';      R = 1/sqrt((1/(1.0*Vn))^2+(obj.a/(0.9*Mn))^2);
                    case 'available';   R = 1/sqrt((1/(Vn/1.50))^2+(obj.a/(Mn/1.67))^2);
                    otherwise; error('Unknown strength type: %s',strength_type)
                end

                strengths    = vertcat(strengths,R);
                limit_states = vertcat(limit_states,'Interaction Strength of the Plate');

                % Lateral-Torsional Buckling of the Plate
                dct = 3;
                Lb  = obj.a;
                Cb  = max([(3+log(Lb/obj.dp))*(1-dct/obj.dp) 1.84]);
                E   = 29000;
                Zp  = obj.tp*obj.dp^2/4;
                Sp  = obj.tp*obj.dp^2/6;
                Mp  = obj.Fy*Zp;
                if Lb*obj.dp/obj.tp^2 <= 0.08*E/obj.Fy
                    Mn  = Mp;
                elseif Lb*obj.dp/obj.tp^2 <= 1.9*E/obj.Fy
                    Mn  = min([Cb*(1.52-0.274*(Lb*obj.dp/obj.tp^2)*obj.Fy/E)*obj.Fy*Sp Mp]);
                else
                    Fcr = 1.9*E*Cb/(Lb*obj.dp/obj.tp^2);
                    Mn  = min([Fcr*Sp Mp]);
                end

                switch lower(strength_type)
                    case 'nominal';     R = Mn/obj.a;
                    case 'design';      R = 0.90*(Mn/obj.a);
                    case 'available';   R = (Mn/obj.a)/1.67;
                    otherwise; error('Unknown strength type: %s',strength_type)
                end

                strengths    = vertcat(strengths,R);
                limit_states = vertcat(limit_states,'Lateral-Torsional Buckling of the Plate');

                % Flexural Rupture of the Plate
                Mn = obj.Fu*obj.Znet;

                switch lower(strength_type)
                    case 'nominal';     R = Mn/obj.a;
                    case 'design';      R = 0.75*(Mn/obj.a);
                    case 'available';   R = (Mn/obj.a)/2.00;
                    otherwise; error('Unknown strength type: %s',strength_type)
                end

                strengths    = vertcat(strengths,R);
                limit_states = vertcat(limit_states,'Flexural Rupture of the Plate');
            end
            
            % Output
            [controlling_R,i] = min(strengths);
                      
            if nargout > 1
                controlling_limit_state = limit_states{i};
            end
            if nargout > 2
                results = struct;
                results.strength_type   = strength_type;
                results.strengths       = strengths;
                results.limit_states    = limit_states;
            end
        end
        function plot_IC_results(obj)
            top_edge    = obj.dp/2;
            bottom_edge = -obj.dp/2;
            right_edge  = (obj.nr-1)*obj.g + obj.leh;
            left_edge   = -obj.a;
            
            % Create BoltGroup objects
            [x,y] = obj.get_bolt_positions; 

            if obj.bolt_hole_def_considered
                Rult = min([obj.Fnv*obj.Ab 2.4*obj.d*obj.tp*obj.Fu]);
                Rult_over_lc = 1.2*obj.tp*obj.Fu;
            else
                Rult = min([obj.Fnv*obj.Ab 3.0*obj.d*obj.tp*obj.Fu]);
                Rult_over_lc = 1.5*obj.tp*obj.Fu;
            end
            
            BG1 = BoltGroup(x,y,Rult);
            
            BG2 = BoltGroupWithTearout(x,y,Rult,Rult_over_lc,obj.dh);
            BG2.edge_top_y       = top_edge;
            BG2.edge_bottom_y    = bottom_edge;
            BG2.edge_right_x     = right_edge;
            
            % Create plot
            fs = figureStyle('Display');
            fs.picture(5,5);
            plot_force_scale = obj.s/Rult;
            
            % Plot bolts and plate boundary
            plot([right_edge right_edge left_edge left_edge right_edge],[top_edge bottom_edge bottom_edge top_edge top_edge],'k')
            scatter(x,y,'k')
            %annotation('arrow',[-1 -1]*obj.e,[0 3])
            ylim([-1 1]*0.6*obj.dp)
            axis equal
            
            % Plot Standard IC results
            [~,IC,rotdir] = BG1.Pn_IC(obj.e,0,0);
            [Rx,Ry] = BG1.Bolt_Forces_IC(IC,rotdir);
            
            my_color1 = 'b';
            scatter(IC(1),IC(2),my_color1)
            for i = 1:length(x)
                plot([IC(1) x(i)],[IC(2) y(i)],'--','Color',my_color1)
                plot(plot_force_scale*[0 Rx(i)]+x(i),plot_force_scale*[0 Ry(i)]+y(i),'-','Color',my_color1)
            end
            
            % Plot Modified IC results
            [~,IC,rotdir] = BG2.Pn_IC(obj.e,0,0);
            [Rx,Ry] = BG2.Bolt_Forces_IC(IC,rotdir);
            
            my_color2 = 'r';
            scatter(IC(1),IC(2),my_color2)
            for i = 1:length(x)
                plot([IC(1) x(i)],[IC(2) y(i)],'--','Color',my_color2)
                plot(plot_force_scale*[0 Rx(i)]+x(i),plot_force_scale*[0 Ry(i)]+y(i),'-','Color',my_color2)
            end            

            
            % % Print out results
            % [Rx,Ry,results] = BG.Bolt_Forces_IC(IC);
            % fprintf('R   angle   lc\n');
            % for i = 1:length(Rx)
            %     fprintf('%.3f   %.3f   %.3f\n',sqrt(Rx(i)^2+Ry(i)^2),rad2deg(results.angle(i)),results.lc(i));
            % end
        end
    end
end


function C = C_from_table(n)
% For eccentricity of: a/2 = 3"/2 = 1.5" when n <= 5
%                      a = 3"            when n >= 6
switch n
    case 2
        C = (1.63+1.18)/2;
    case 3
        C = (2.71+2.23)/2;
    case 4
        C = (3.75+3.32)/2;
    case 5
        C = (4.77+4.39)/2;
    case 6
        C = 4.98;
    case 7
        C = 6.06;
    case 8
        C = 7.12;
    case 9
        C = 8.17;
    case 10
        C = 9.21;
    case 11
        C = 10.2;
    case 12
        C = 11.3;
    otherwise
        error('Value not stored');
end
end
