function [V,qR,n,re,lam]=fall_speed_V_qR_iteration(P,V)
%P in mm/hr - this stays constant
%V is the best guess for V
%Iterates to solve for V and qR

c = 523.6;
rhoa=1;

diff = 1e9;
tol = 0.01e-3; %qR tolerance in kg/kg
n=0;
qR = 1e9;

while max(diff)>tol & n<20
    qR_old = qR;
    n=n+1;

    qR=P/3600 /rhoa ./ V;  %R = rhoa * qR * V
    [V,lam,re] = fall_speed_mass_weighted(qR,P);

    diff = abs(qR - qR_old);
end

