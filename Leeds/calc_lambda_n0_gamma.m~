function [n0,lam] = calc_lambda_n0_gamma(m1, m2, p1, p2, mu)
%Based on subroutine get_lam_n0_2M(m1, m2, p1, p2, mu, lam, n0) in
%lookup.F90 from CASIM code

%    real(wp), intent(in) :: m1, m2, mu
%    real(wp), intent(in) :: p1, p2
%    real(wp), intent(out) :: lam, n0

%    real(wp) :: p, m
%    real(wp) :: j1

    j1=1.0/(p1-p2);

    lam=((Gamma(1.0+mu+p1)./Gamma(1.0+mu+p2)).*(m2./m1)).^(j1); %Note here that the units of m1 and m2 cancel, so whether they are in per m3 or per kg does not matter for the lambda

    m=m2;
    p=p2;

    %N.B.- the n0 here is just the number concentration - this is how it is
    %coded in the UM... 
    %Need to use the size dist in the form :-
    % ND = n0.* ( lam.^(1+mu) / gamma(1+mu) ) .*D.^mu.*exp(-lam.*D);
    n0=lam.^(p).*m.*Gamma(1.0+mu)./Gamma(1.0+mu+p);

