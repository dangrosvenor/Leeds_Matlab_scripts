function [h]=m_dpcolor(X,Y,C)
%function [h]=m_dpcolor(X,Y,C)
%a wrapper for pcolor that allows fields to be plotted where the size of C is one row
%and column less than that of X and Y, so that the last row and column of C
%are actually displayed (Matlab discards them and there doesn't seem to be
%a way to make them get plotted, even though it is stated in the help doc
%that it should plot in this situation)
%So, this just adds dummy rows and columns to C to make it the same size as
%X and Y

C(:,end+1)=NaN*ones([size(C,1) 1]);
C(end+1,:)=NaN*ones([1 size(C,2)]);

h=m_pcolor(X,Y,C);
