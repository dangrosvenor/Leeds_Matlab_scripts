function [me,std] = pdf_me_std(x,F)
%function [me,std] = pdf_me_std(x,F)
%calculates the mean and standard deviation of a given 1D PDF (can be
%normalised or not) described by mid-points x and frequencies F

if size(x,1)==1
    x=x';
end

if size(F,1)==1
    F=F';
end


N = sum(F);
me = sum(F.*x)./N;
std = sqrt (   sum( F.*(x-me).^2 , 1 ) ./ (N)   );