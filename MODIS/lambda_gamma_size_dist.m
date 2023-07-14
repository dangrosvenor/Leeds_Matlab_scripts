function [lam,n0] = lambda_gamma_size_dist(qR,alpha,na,nb,c,d,rhoa,nR)
%function lam = lambda_gamma_size_dist(qR)
%qR in kg/kg 
%Returns lamb=da in units m^-1

n0=-999;

dist='UM CASIM 2M';

switch dist
    case 'UM CASIM 2M'
        %Based on get_lam_n0_2M in lookup.F90
        m1 = qR; %mass
        m2 = nR; %number
        mu = alpha;
        
        p1=3; p2=0; p3=6; %from rain_params bit of code
        j1=1./(p1-p2);

        lam = ((gamma(1.+mu+p1)./gamma(1.+mu+p2)) .*(m2./m1)).^(j1);
       
        m=m2;
        p=p2;

        n0 = lam.^p.*m.*gamma(1.+mu)./gamma(1.+mu+p);
        
    otherwise

if ~exist('alpha')
    dist = 'gamma';
    dist = 'Marshall';

    %constants for rain

    nb = 0;  %This makes n(D) of the form na * D^alpha * exp(-lambda*D)
    c = 523.6;
    d = 3;
    rhoa = 1;

    switch dist
        case 'gamma'
            alpha = 2.5;
            na = 1.1e15;
        case 'Marshall'
            alpha = 0;
            na = 8e6;
    end

end


lam = ( na*c*gamma(1+alpha+d) ./ (rhoa*qR) ) .^( 1/(1+alpha+d-nb) ) ;

end