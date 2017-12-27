classdef BoltGroup
    
    properties
        x
        y
        Rult
        Dmax = 0.34;
    end
    
    methods
        function obj = BoltGroup(x,y,Rult)
            if nargin < 3
                Rult = 1;
            end
            assert(isequal(size(x),size(y)),'x and y should be the same size')
            obj.x = x;
            obj.y = y;
            obj.Rult = Rult;
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
        function [Rx,Ry] = Bolt_Forces_IC(obj,IC)
            ICx = IC(1);
            ICy = IC(2);
            
            % Calculate r values
            rx = obj.x-ICx;
            ry = obj.y-ICy;
            r = sqrt(rx.^2 + ry.^2);
            
            % Compute delta
            delta = obj.Dmax*(r/max(r));
            
            % Compute strengths of bolts
            R = obj.Rult*(1-exp(-10*delta)).^0.55;
            
            % Component breakdown
            angle = atan2(ry,rx) + pi/2;
            Rx = cos(angle).*R;
            Ry = sin(angle).*R;
        end        
        
        function [Pn,IC] = Pn_IC(obj,xP,yP,theta)
            % theta input as degrees clockwise from y-axis
            
            % convert theta to radians counter-clockwise from x-axis
            theta = deg2rad(270-theta);
            
            options = struct;
            options.Algorithm = 'levenberg-marquardt';
            % options.Display = 'iter-detailed';
            options.Display = 'off';
            
            [CGx,CGy] = obj.centroid();
            ICo = [CGx CGy];
            [IC,~,exitflag,output] = fsolve(@(IC)error_Pn_IC(obj,xP,yP,theta,IC),ICo,options);
            if exitflag <= 0
                fprintf(output);
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
            options.Display = 'off';
            
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
                    plot([0 Rx(i)]+obj.x(i),[0 Ry(i)]+obj.y(i),'-r')
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
error = sqrt(errorx^2 + errory^2);

end
