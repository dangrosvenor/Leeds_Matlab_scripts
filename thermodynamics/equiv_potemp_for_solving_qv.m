function y=equiv_potemp_for_solving_qv(qv,TK,P,equiv_potemp_val)
%function y=equiv_potemp_for_solving_qv(qv,TK,P,equiv_potemp_val)
%function y to be solved for y=0
%y=equivalent_potemp(TK,P,qv)-equiv_potemp_val;
%Given a qv value to be varied, as well as fixed TK and P values
%fzero can be used to find the qv value that is consistant with 
%the supplied equiv_potemp_val value.
%e.g. fzero(@equiv_potemp_for_solving_qv,[0 5e-3],[],280,800e2,304.18)
%returns qv=0.002 since equivalent_potemp(TK=280,P=800e2,qv)=304.18
%[0 5e-3] is the range within which the answer lies - important to include
%this as otherwise a value might not be found
y=equivalent_potemp(TK,P,qv)-equiv_potemp_val;