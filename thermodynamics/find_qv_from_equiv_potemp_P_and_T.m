%finds the vapour mixing ratio given the equiv potemp, pressure and temperature using fzero
function [qv]=find_qv_from_equiv_potemp_P_and_T(equiv,P,TK)

for i=1:length(equiv)
    qv(i)=fzero(@equiv_potemp_for_solving_qv,[0 4e-3],[],TK(i),P(i),equiv(i));
end