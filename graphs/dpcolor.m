function [h,hback]=dpcolor(X,Y,C)
%function [h]=dpcolor(X,Y,C)
%a wrapper for pcolor that allows fields to be plotted where the size of C is one row
%and column less than that of X and Y, so that the last row and column of C
%are actually displayed (Matlab discards them and there doesn't seem to be
%a way to make them get plotted, even though it is stated in the help doc
%that it should plot in this situation)
%So, this just adds dummy rows and columns to C to make it the same size as
%X and Y

if nargin==1
    C=X;
end

C(:,end+1)=NaN*ones([size(C,1) 1]);
C(end+1,:)=NaN*ones([1 size(C,2)]);

minC = minALL(C);
maxC = maxALL(C);

if nargin==1
    hback2=pcolor(ones(size(C))); %if add two backgrounds on then it reduces the aliasing hatching
    %in the final PDF
    hold on
    hback=pcolor(ones(size(C)));         %
    hold on
    h=pcolor(C);
%    shading flat
    set(hback,'facecolor',[0 0 0]);
else
    hback2=pcolor(X,Y,ones(size(C))); %if add two backgrounds on then it reduces the aliasing hatching
    %in the final PDF
    hold on
    hback=pcolor(X,Y,ones(size(C)));   
%    shading flat
    hold on
    h=pcolor(X,Y,C);
%    shading flat
    set(hback,'facecolor',[0 0 0]);
end

if maxC>minC
    caxis([minC maxC]);
end