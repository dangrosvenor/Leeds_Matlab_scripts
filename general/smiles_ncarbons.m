function [ncarbons]=smiles_ncarbons(str_out)
%counts the number of carbons in the base chain (not including oxygen atoms etc.)

nleft=9999;
while nleft~=0
    [str_out,str_in,nleft,nleft2,npos,npos_start]=smiles_bracket(str_out); 
    %basically keeps removing complete brackets until are left with a chain of standalone carbons and other non-bracketed atoms (e.g. 'COCO')
end

%then count the number of carbons left
 posC=findstr('C',str_out); 
 ncarbons=length(posC);

