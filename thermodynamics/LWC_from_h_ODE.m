function L = LWC_from_h_ODE(h,T_cb,P_cb,power)
%returns L in g/m3 for an adiabatic ascent from cloud base temp of T_cb and
%pressure = P_cb over a distance of h (cloud thickness = h)

if nargin<4
    power=1;
end

for i=1:length(h)
    if h(i)==0
        L(i)=0;
    else
        [z,y]=ode45(@ode_test_multi_dLdz,[0 h(i)],[0 T_cb P_cb]);
        %E.g. %  [z,y]=ode45(@ode_test_multi_dLdz,[0 h],[0 283 950e2]);
        % Here 0 = initial L = 0, 283 = T_cloudbase and 950e2 = P_cloudbase
        % z is a vector of the z values that correspond to the values returned in y
        % (size = [M 1])
        % y is an array of size [M 3] where y(:,i) are the solution values for
        % variable i

        L(i) = (y(end,1)).^power;  %kg/m3

    end

end