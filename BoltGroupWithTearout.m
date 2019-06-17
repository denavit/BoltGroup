classdef BoltGroupWithTearout
    
    properties
        x
        y
        Rult
        Rult_over_lc
        dh
        Dmax = 0.34;
        
        edge_right_x    = Inf;
        edge_left_x     = -Inf;
        edge_top_y      = Inf;
        edge_bottom_y   = -Inf;
        
        plot_force_scale = 1;
        ignore_clear_distance_between_holes = false;
        
        fsolve_Display = 'off';        
    end
    
    methods
        function obj = BoltGroupWithTearout(x,y,Rult,Rult_over_lc,dh)
            assert(isequal(size(x),size(y)),'x and y must be arrays of the same size')
            assert(isscalar(Rult)==1,'Rult must be a scalar');
            assert(isscalar(Rult_over_lc)==1,'Rult_over_lc must be a scalar');
            assert(isscalar(dh)==1,'dh must be a scalar');
            obj.x = x;
            obj.y = y;
            obj.Rult = Rult;
            obj.Rult_over_lc = Rult_over_lc;
            obj.dh = dh;
        end
        
        function [xc,yc] = centroid(obj)
            xc = mean(obj.x);
            yc = mean(obj.y);
        end
        function [Rx,Ry,results] = Bolt_Forces_IC(obj,IC,d)
            % This function returns the force each of the bolts imparts 
            % on the plate.
            %
            % Direction is the sign of the rotation of the plate with
            % respect to the bolt group (counterclockwise positive).
            %
            
            ICx = IC(1);
            ICy = IC(2);
            
            % Check edges
            assert(all(obj.x+obj.dh/2<obj.edge_right_x) ,'All bolts must be left of right edge');
            assert(all(obj.x-obj.dh/2>obj.edge_left_x)  ,'All bolts must be right of left edge');
            assert(all(obj.y-obj.dh/2>obj.edge_bottom_y),'All bolts must be above bottom edge');
            assert(all(obj.y+obj.dh/2<obj.edge_top_y)   ,'All bolts must be below top edge');
            
            % Calculate r values
            rx = obj.x-ICx;
            ry = obj.y-ICy;
            r = sqrt(rx.^2 + ry.^2);
            
            % Compute delta
            delta = obj.Dmax*(r/max(r));
            
            % Compute clear distance
            angle = atan2(ry,rx) - sign(d)*pi/2;
            lc = inf(size(obj.x));
            for i = 1:length(obj.x)
                % Check against right edge
                if (obj.edge_right_x - obj.x(i)) > 0 && cos(angle(i)) > 0
                    lci = (obj.edge_right_x - obj.x(i))/cos(angle(i)) - obj.dh/2;
                    if lci < lc(i)
                        lc(i) = lci;
                    end
                end
                
                % Check against left edge
                if (obj.edge_left_x - obj.x(i)) < 0 && cos(angle(i)) < 0
                    lci = (obj.edge_left_x - obj.x(i))/cos(angle(i)) - obj.dh/2;
                    if lci < lc(i)
                        lc(i) = lci;
                    end
                end
                
                % Check against top edge
                if (obj.edge_top_y - obj.y(i)) > 0 && sin(angle(i)) > 0
                    lci = (obj.edge_top_y - obj.y(i))/sin(angle(i)) - obj.dh/2;
                    if lci < lc(i)
                        lc(i) = lci;
                    end
                end
                
                % Check against bottom edge
                if (obj.edge_bottom_y - obj.y(i)) < 0 && sin(angle(i)) < 0
                    lci = (obj.edge_bottom_y - obj.y(i))/sin(angle(i)) - obj.dh/2;
                    if lci < lc(i)
                        lc(i) = lci;
                    end
                end
                
                % Check against other bolts
                if ~obj.ignore_clear_distance_between_holes
                    for j = 1:length(obj.x)
                        if i ~= j
                            lci = clear_distance_between_bolt_holes(obj.x(i),obj.y(i),obj.x(j),obj.y(j),obj.dh,angle(i));
                            if lci < lc(i)
                                lc(i) = lci;
                            end
                        end
                    end
                end
            end
            
            % Compute strengths of bolts
            Rult_each = min(obj.Rult,lc*obj.Rult_over_lc);
            R = Rult_each.*(1-exp(-10*delta)).^0.55;
            
            % Component breakdown
            Rx = cos(angle).*R;
            Ry = sin(angle).*R;
            
            % Extra output
            if nargout > 2
                results = struct;
                results.angle   = angle;
                results.lc      = lc;
            end
        end        
        
        function [Pn,IC,d] = Pn_IC(obj,xP,yP,theta)
            % theta input as degrees clockwise from y-axis
            
            % convert theta to radians counter-clockwise from x-axis
            theta = deg2rad(270-theta);
            
            % Find instantaneous center
            options = struct;
            options.Display = obj.fsolve_Display;
            
            [CGx,CGy] = obj.centroid();
            ICo = [CGx CGy];
            [IC,~,exitflag,output] = fsolve(@(IC)error_Pn_IC(obj,IC,xP,yP,theta),ICo,options);
            if exitflag <= 0
                ICo = [CGx-3 CGy];
                [IC,~,exitflag,output] = fsolve(@(IC)error_Pn_IC(obj,IC,xP,yP,theta),ICo,options);
                if exitflag <= 0
                    fprintf('%s\n',output.message);
                    error('Unable to determine the instantaneous center');
                end
            end
            
            % Compute Pn from instantaneous center 
            d = direction(IC,xP,yP,theta);
            [Rx,Ry] = obj.Bolt_Forces_IC(IC,d);
            Pn = sqrt(sum(Rx)^2 + sum(Ry)^2);
            
            % Output
            if nargout < 2
                clear IC
            end
            if nargout < 3
                clear d
            end
        end
        function [Mn,IC] = Mn_IC(obj,d)
            
            % Find instantaneous center
            options = struct;
            options.Display = obj.fsolve_Display;
            
            [CGx,CGy] = obj.centroid();
            ICo = [CGx CGy];
            [IC,~,exitflag,output] = fsolve(@(IC)error_Mn_IC(obj,IC,d),ICo,options);
            if exitflag <= 0
                fprintf(output);
                error('Unable to determine the instantaneous center');
            end
            
            % Compute Mn from instantaneous center
            [Rx,Ry] = obj.Bolt_Forces_IC(IC,d);
            Mn = abs(sum(-Rx.*obj.y) + sum(Ry.*obj.x));
            
            % Output
            if nargout < 2
                clear IC
            end            
        end        
        function plot(obj,IC,xP,yP,theta)
            hold all
            scatter(obj.x,obj.y,'k')
            
            [xc,yc] = obj.centroid;
            scatter(xc,yc,'r')
            
            if nargin > 1
                if nargin == 3
                    d = xP;
                else
                    d = direction(IC,xP,yP,theta);
                end                
                scatter(IC(1),IC(2),'b')
                [Rx,Ry] = obj.Bolt_Forces_IC(IC,d);
                for i = 1:length(obj.x)
                    plot([IC(1) obj.x(i)],[IC(2) obj.y(i)],'--k')
                    plot(obj.plot_force_scale*[0 Rx(i)]+obj.x(i),obj.plot_force_scale*[0 Ry(i)]+obj.y(i),'-r')
                end
            end
        end
        
    end
end

function d = direction(IC,xP,yP,theta)
d = sign((xP-IC(1))*sin(theta) - (yP-IC(2))*cos(theta));
end

function error = error_Pn_IC(obj,IC,xP,yP,theta)   
ICx = IC(1);
ICy = IC(2);
d = direction(IC,xP,yP,theta);
[Rx,Ry] = obj.Bolt_Forces_IC(IC,d);

% Calculate r values
rx = obj.x-ICx;
ry = obj.y-ICy;

% Sum Moments to solve for P
P = sum(ry.*Rx - rx.*Ry)/((xP-IC(1))*sin(theta) - (yP-IC(2))*cos(theta));

% Sum forces
errorx = sum(Rx)+P*cos(theta);
errory = sum(Ry)+P*sin(theta);
error = [errorx errory];
end

function error = error_Mn_IC(obj,IC,d)
[Rx,Ry] = obj.Bolt_Forces_IC(IC,d);
error = [sum(Rx) sum(Ry)];
end