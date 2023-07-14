function[SSw, SSi, nTa, x, bound, line_data]...
    =  CFDFlow(P, Tw, Tc, Qtot, Qs, Ta)

% CFDFLOW Calculate lamina conditions such as the temperature and
% supersaturations.
%
%   [SSw, SSi, nTa, x, bound, line_data] =  CFDFlow(P, Tw, Tc, Qtot, Qs, Ta)
%
% CFDFlow is a function that calculates variables associated with the
% lamina position as a function of the measured CFDC variables.  The inputs
% are;
%   (1) P = the chamber pressure, mb
%   (2) Tw = the average warm wall temperature, degC
%   (3) Tc = the average cold wall temperature, degC
%   (4) Qtot = the total flow rate, lpm
%   (5) Qs = the sheath flow rate, lpm
% Optional inputs
%   (6) Ta = previously calculated aerosol lamina temperature, degC; used
%           to convert mass flows to volumetric flows.  If this argument is
%           not specified, the flows are assumed to be volumetric
%
% The function outputs are:
%   (1) SSw = the calculated supersaturation wrt water, %
%   (2) SSi = the calculated supersaturation wrt ice, %
%   (3) nTa = the new lamina temperature, degC
%   (4) x = the lamina position; this is the ratio of the distance of the
%           lamina from the cold wall to the annular gap width, no units
%   (5) bound = the conditions at the boundaries of the lamina. 'bound' is
%               a 3x2 matrix with the first column representing the SSw,
%               SSi and T at the boundary closest to the cold wall and the
%               second column representing the values closest to the warm
%               wall.
%
% ALL UNITS ARE CGS
%
%
% See also ESATW, ESATI, FZERO.

%   Copyright 2007 Colorado State University
%   $Revision: 2.0 $  $Date: 2007/08/03 $

%% Check the number of arguments and initialize the necessary values

% Only five arguments are required - if they are not there, throw an error
disp('Hello')

if nargin < 5
    error('MATLAB:CFDFlow:NotEnoughInputs',...
        'Not enough input arguments.  See CFDFLOW.');
elseif nargin > 6
    error('MATLAB:CFDFlow:TooManyInputs',...
        'Too many input arguments.  See CFDFLOW.');
end
disp('Hello2')
% Change Qs to the sample flow
Qs = Qtot - Qs;

if nargin == 6
    % Conversion constant for going from slpm -> lpm
    C = 1013.25/P*(Ta+273.15)/273.15;

    % Convert the mass flow rates to volumetric flow rates given
    Qtot = Qtot*C;
    Qs   = Qs*C;
end

% Gap half-width in cm
d = 0.5*1.12;
WidthMean = 28.97;
disp('Hello3')
% Square and cube of gap half-width; these are used throughout, so we do
% them once here
d2 = d^2;
d3 = d^3;

% Dry air gas constant
Rd = 2.87e6;

% Gravitational acceleration
g = 980.6;

% Convert the pressure to dyne/decimeter^2 from mb
P = P*1000;

% Fraction of flow that is sample
f = Qs/Qtot;

% Assume that the flow is evenly divided - fas is the cold side fraction
% and fbs is the warm side fraction
fas = 0.5*(1-f);
fbs = 1 -f -fas;
disp('Hello4')
% Calculate the dynamic viscocity
mu = 1.0e-4*(1.718+(0.005*(0.5*(Tw+Tc))));

Gam = P*g*d2/(12.0*mu*Rd*(273.16+(0.5*(Tw+Tc)))^2);

% The temperature difference between the warm and cold walls
delT = Tw - Tc;
disp('Hello5')
% Convert the flow rate from lpm -> cc/s
Q = Qtot*1000/60;
disp('Hello6')
% This is the mean velocity in cm/s
velocmean = Q/(2*d*WidthMean);
disp('Hello7')
% Critical flow rate in cm/s
Qcrit = 2*d*WidthMean/3*Gam*delT;

if (Qcrit > Q)
    warning('MATLAB:CFDFLOW:CriticalFlow',...
        'Total flow is less than critical flow. See CFDFlow.');
end

% This is the range in which to find zeros
r = [-d,d];
disp('Hello8')
%% Find the lamina position

% This is the a value (closest to cold wall)
x(2,1) = fzero(@fx, r, [], d, P);
disp('Hello9')
% This is the b value (closest to warm wall)
x(3,1) = fzero(@gx, r);

if x(3) < x(2)
    warning('MATLAB:CFDFLOW:BLTA', 'b < a. See CFDFlow.');
end

% This is the half-width of the annular gap
x(4,1) = d;

% Calculate x, average fractional position of the aerosol lamina.
% This fractional position is relative to the cold wall.
x(1,1) = (x(2)+x(3)+2*d)/(4*d);

%% Calculate the average values at the lamina center
nTa = Tc + x(1)*delT;

ec = esati(Tc, 'Murph');
ew = esati(Tw,'Murph');

e = ec + x(1)*(ew - ec);

SSw = e/esatw(nTa, 'Murph')*100 - 100;

SSi = e/esati(nTa,'Murph')*100 - 100;

% Information at the boundaries
ea = ec+ (x(2)+d)/(2*d)*(ew-ec);
eb = ec+ (x(3)+d)/(2*d)*(ew-ec);
aTa = Tc + (x(2)+d)/(2*d)*delT;
bTa = Tc + (x(3)+d)/(2*d)*delT;
SSwa = ea/esatw(aTa, 'Murph')*100 - 100;
SSwb = eb/esatw(bTa, 'Murph')*100 - 100;
SSia = ea/esati(aTa,'Murph')*100 - 100;
SSib = eb/esati(bTa,'Murph')*100 - 100;

bound = [SSwa SSwb;SSia SSib; aTa bTa];

z = [0:0.05:1]';
Tz = Tc + z*delT;
ez = ec + z*(ew - ec);
SSwz = ez./esatw(Tz, 'Murph')*100 - 100;

SSiz = ez./esati(Tz,'Murph')*100 - 100;

line_data = [z*d*2 Tz SSwz SSiz];


%% These are functions that are used to find the boundaries of the lamina
    %function flux = Flux(p,velocmean)
    function flux = Flux(p,velocmean)
        z = p^2;
disp('Hello10')
        flux = 1.5*velocmean*(p-p*z/(3*d^2)) + Gam*delT*(z^2/(4*d3)-z/(2*d));
        disp('Hello11')
    end

    % The following functions are used to find the boundaries of the lamina
    % - they should be zeros at their respective boundaries.
    function F = fx(a,d,p)
        
        F = Flux(p,a) - Flux(p,-d) - fas*Q/WidthMean;
    end

    function G = gx(b)
        
        G = Flux(d) - Flux(b) - fbs*Q/WidthMean;
        
    end

end