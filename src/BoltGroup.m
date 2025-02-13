classdef BoltGroup
    
    properties
        x
        y
        Rult
        Dmax = 0.34;
        load_deformation_type = 'standard';
        Rslip = 0;
        Dslip = 0;
        
        plot_force_scale = 1;
        
        fsolve_Display = 'off';
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
        function [Rx,Ry,results] = Bolt_Forces_IC(obj,IC,d)
            % This function returns the force each of the bolts imparts 
            % on the plate.
            %
            % Direction is the sign of the rotation of the plate with
            % respect to the bolt group (counterclockwise positive).
            %
            
            ICx = IC(1);
            ICy = IC(2);
            
            % Calculate r values
            rx = obj.x-ICx;
            ry = obj.y-ICy;
            r = sqrt(rx.^2 + ry.^2);
            
            % Compute force in bolts
            switch lower(obj.load_deformation_type)
                case 'standard'
                    delta = obj.Dmax*(r/max(r));
                    R = obj.Rult*(1-exp(-10*delta)).^0.55;
                case 'standard_with_slip'
                    delta = (obj.Dmax+obj.Dslip)*(r/max(r));
                    Dslip0 = -(1/10)*log(1-(obj.Rslip/obj.Rult)^(1/0.55));
                    
                    R = nan(size(delta));
                    for i = 1:length(delta)
                        if delta(i) <= Dslip0
                            R(i) = obj.Rult*(1-exp(-10*delta(i))).^0.55;
                        elseif delta(i) <= Dslip0 + obj.Dslip
                            R(i) = obj.Rslip;
                        else
                            R(i) = obj.Rult*(1-exp(-10*(delta(i)-obj.Dslip))).^0.55;
                        end
                    end
                case 'elastic'
                    R = obj.Rult*(r/max(r));
                case 'plastic'
                    %R = obj.Rult;
                    R = min(100*obj.Rult*(r/max(r)),obj.Rult);
                otherwise
                    error('Unknown load deformation type: %s',obj.load_deformation_type);
            end
            
            % Component breakdown
            angle = atan2(ry,rx) - sign(d)*pi/2;
            Rx = cos(angle).*R;
            Ry = sin(angle).*R;
            
            % Extra output
            if nargout > 2
                results = struct;
                results.angle   = angle;
                results.lc      = nan(size(angle));
            end
        end        
        function [Pn,IC,d] = Pn_IC(obj,xP,yP,theta)
            % theta input as degrees clockwise from y-axis
            
            % Convert theta to radians counter-clockwise from x-axis
            theta = deg2rad(270-theta);
            
            % Find instantaneous center
            options = struct;
            options.Display = obj.fsolve_Display;
            
            [CGx,CGy] = obj.centroid();
            ICo = [CGx CGy];
            [IC,~,exitflag,output] = fsolve(@(IC)error_Pn_IC(obj,IC,xP,yP,theta),ICo,options);
            if exitflag <= 0
                fprintf('%s\n',output.message);
                error('Unable to determine the instantaneous center');
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
            if nargin < 3
                d = 1;
            end
            
            % Find instantaneous center
            options = struct;
            options.Display = obj.fsolve_Display;
            
            [CGx,CGy] = obj.centroid();
            ICo = [CGx CGy];
            [IC,~,exitflag,output] = fsolve(@(IC)error_Mn_IC(obj,IC,d),ICo,options);
            if exitflag <= 0
                fprintf('%s\n',output.message);
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