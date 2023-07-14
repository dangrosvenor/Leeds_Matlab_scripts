function [pdfarr] = ukca_ACTIVATE_PDF_FUNC(array,sig,me,alpha)

for k=1:length(array)
    
    tmp1(k)= -((array(k)-me).^2./(2.*sig.^2));
    tmp1_out(k) = exp(tmp1(k));

    pdfarr(k)=(1.0./((2.0.*pi).^0.5)).*(1.0./sig).*                    ...
    tmp1_out(k) .*                                                    ...
    (1+erf(alpha.*(array(k)-me).*                        ...
    (1.0./((2.0).^0.5)).*(1.0./sig)));
          
end