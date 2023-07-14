function [t] = signifcance_of_diff_two_PDFs()
%

Spooled = sqrt ( (s1.^2.*(n1-1) + s2.^2.*(n2-1)) ./ (n1+n2-2) );

t = (me01 - me02)./Spooled .* sqrt( (n1.*n2) ./ (n1+n2) );