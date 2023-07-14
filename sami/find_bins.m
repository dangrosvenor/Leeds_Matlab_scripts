function [n_bins,n_gases,rlen]=find_bins(vec)
% Funktio etsii binien ja kaasujen lukumäärät, kun data on muodossa
%   vec=[t (s), [r_bin(1:n_bins) (m)], RH, T(K), [c_gas(1:n_gases) (mol/cm^3)]]
% rlen on yhden ajanhetken tietojen pituus
%
% vec(t)=[time,r_bin(n_bin),RH,T,c_gas(n_gas)]
for i=2:length(vec)
    aa=str2num(vec{i});
    if (aa>0.1) % RH
        n_bins=i-2
        break
    end
end
%
for i=0:length(vec)
    aa=str2num(vec{n_bins+4+i});
    if (aa>0.01)
        n_gases=i
        break
    end
end
rlen=n_gases+n_bins+3;
%--------------------------------------------------------------------------