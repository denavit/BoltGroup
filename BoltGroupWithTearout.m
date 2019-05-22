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
        function P = P_IC(obj,IC)   
            [Rx,Ry] = obj.Bolt_Forces_IC(IC);
            P = sqrt(sum(Rx)^2 + sum(Ry)^2);
        end
        function M = M_IC(obj,IC)
            [Rx,Ry] = obj.Bolt_Forces_IC(IC);
            M = sum(-Rx.*obj.y) + sum(Ry.*obj.x);
        end
        function [Rx,Ry,results] = Bolt_Forces_IC(obj,IC)
            ICx = IC(1);
            ICy = IC(2);
            
            % Calculate r values
            rx = obj.x-ICx;
            ry = obj.y-ICy;
            r = sqrt(rx.^2 + ry.^2);
            
            % Compute delta
            delta = obj.Dmax*(r/max(r));
            
            % Compute clear distance
            angle = atan2(ry,rx) + pi/2;
            rad2deg(angle);
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
        
        function [Pn,IC] = Pn_IC(obj,xP,yP,theta)
            % theta input as degrees clockwise from y-axis
            
            % convert theta to radians counter-clockwise from x-axis
            theta = deg2rad(270-theta);
            
            options = struct;
            options.Display = obj.fsolve_Display;
            
            [CGx,CGy] = obj.centroid();
            ICo = [CGx CGy];
            [IC,~,exitflag,output] = fsolve(@(IC)error_Pn_IC(obj,xP,yP,theta,IC),ICo,options);
            if exitflag <= 0
                fprintf('%s\n',output.message);
                error('Unable to determine the instantaneous center');
            end
            Pn = obj.P_IC(IC);
            if nargout < 2
                clear IC
            end
        end
        function [Mn,IC] = Mn_IC(obj)
            options = struct;
            options.Algorithm = 'levenberg-marquardt';
            options.Display = obj.fsolve_Display;
            
            [CGx,CGy] = obj.centroid();
            ICo = [CGx CGy];
            [IC,~,exitflag,output] = fsolve(@(IC)obj.P_IC(IC),ICo,options);
            if exitflag <= 0
                fprintf(output);
                error('Unable to determine the instantaneous center');
            end
            Mn = abs(obj.M_IC(IC));
            if nargout < 2
                clear IC
            end
        end        
        

        
        function plot(obj,IC)
            hold all
            scatter(obj.x,obj.y,'k')
            
            [xc,yc] = obj.centroid;
            scatter(xc,yc,'r')
            
            if nargin > 1
                scatter(IC(1),IC(2),'b')
                [Rx,Ry] = obj.Bolt_Forces_IC(IC);
                for i = 1:length(obj.x)
                    plot([IC(1) obj.x(i)],[IC(2) obj.y(i)],'--k')
                    plot(obj.plot_force_scale*[0 Rx(i)]+obj.x(i),obj.plot_force_scale*[0 Ry(i)]+obj.y(i),'-r')
                end
            end
        end
        
    end
end




function error = error_Pn_IC(obj,xP,yP,theta,IC)   

ICx = IC(1);
ICy = IC(2);
[Rx,Ry] = obj.Bolt_Forces_IC(IC);
R = sqrt(Rx.^2+Ry.^2);

% Calculate r values
rx = obj.x-ICx;
ry = obj.y-ICy;
r = sqrt(rx.^2 + ry.^2);

% Sum Moments to solve for P
a = sin(theta);
b = -cos(theta);
c = cos(theta)*yP-sin(theta)*xP;
e = (a*ICx + b*ICy + c)/sqrt(a^2+b^2); % @todo - should this be pos, neg, or abs value????
P = sum(r.*R)/e;

% Sum forces
errorx = sum(Rx)+P*cos(theta);
errory = sum(Ry)+P*sin(theta);
error = [errorx errory];

end
