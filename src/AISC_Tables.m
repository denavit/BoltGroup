function [C,Cp] = AISC_Tables(table,angle,s,ex,n)

% Define Bolt Geometry
xi = ones(1,n);
yi = 0:s:((n-1)*s);
switch table
    case '7-6'
        x = 0*xi;
        y = yi;
        
    case '7-7'
        x = [0*xi 3*xi];
        y = [yi yi];

    case '7-8'
        x = [0*xi 5.5*xi];
        y = [yi yi];

    case '7-9'
        x = [0*xi 8*xi];
        y = [yi yi];
       
    case '7-10'
        x = [0*xi 3*xi 6*xi];
        y = [yi yi yi];
       
    case '7-11'
        x = [0*xi 6*xi 12*xi];
        y = [yi yi yi];

    case '7-12'
        x = [0*xi 3*xi 6*xi 9*xi];
        y = [yi yi yi yi];
        
    case '7-13'
        x = [0*xi 4*xi 8*xi 12*xi];
        y = [yi yi yi yi];
        
    otherwise
        error('Unknown table: %s',table)
end

% Define Bolt Group Object
BG = BoltGroup(x-mean(x),y-mean(y));

% Solve for the strength coefficient
C = BG.Pn_IC(ex,0,angle);

% Solve for the moment coefficient (if requested)
if nargout > 1
    Cp = BG.Mn_IC();
end

end

