%estimate effect of using cw as a function of temperature all throughout
%the cloud

%From MODIS we have the cloud top (z=H, z=0 is cloud base) re:-
% re(H) = 3*L(H)./(rhoW.*A)      --- (1)    [   A = 4*pi*(N.*k).^(1/3) .* (3.*L(H)./(4*pi*rhoW)).^(2/3);   ]
%           N is assumed constant with height
% Procedure is to guess H (cloud thickness). Then can calculate an N value
% using (1)
% Once have N and H can calculate what the integral I in the tau equation
% needs to be :-
% tau = pi.*Q.*(N.*k).^(1/3).*(3./(4*pi*rhoW)).^(2/3) .* I    --- (2)
%       I = integral{0,H} ( L(h).^(2/3) ) dh  --- (3)
%       where L(h) = integral{0,h} ( cw(h*) ) dh*   -- cw is the
%       condensation rate, which we know as a function of (T,P) and
%       therefore h (from moist adiabatic lapse rate).
% We know I and so can probably estimate the H value required to produce
% it.
% Then can use the new H to reiterate.

% Base the inital guess of H on the cloud top cw

h = 


