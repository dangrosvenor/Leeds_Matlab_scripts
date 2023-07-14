function y = mmfsolve_test_Dan(x)
%x is a vector with values for two different variables (e.g. x(1)=a,
%x(2)=b)

% 
% y(1,:) = x(1,:).^2 + 2.*x(2,:);
% %y(2,:) = 2.*x(1,:) + x(2,:).^3 - 100;
% y(2,:) = 0; %just a plane surface at z=0

%Think doesn't need to be a matrix (one value called at a time for each variable)
y(1) = x(1).^2 + 2.*x(2);
%y(2) = 2.*x(1) + x(2).^3 - 100;
y(2) = 0; %just a plane surface at z=0

y=y';

