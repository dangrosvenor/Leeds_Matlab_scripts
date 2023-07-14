err_str = {'cw','fad','tau','k','re','Vert strat (Odran)'};

err_sq = [0.25*[24 30 36].^2 16^2 25/4*28^2 20^2]
total = sqrt(sum(err_sq))
contrib = 100* err_sq./sum(err_sq)


err_sq_nore = [0.25*[24 30 36].^2 16^2 25/4*0^2 20^2];
total_nore = sqrt(sum(err_sq_nore))
contrib_nore = 100 * err_sq_nore./sum(err_sq_nore)

fad=0.8; dfad=0.1; pdfad = 100*dfad/fad;
pdfad=30; %From Merk - mean of 0.66 std dev of 0.22 giving approx 30% (actually 33%)
%k uncertainty of 12.5% from Merk too
%17% chosen for re based on VOCALS and LES studies.
%What to choose for tau? Perhaps 17% as for re for now??

%err_sq_Dan = [0.25*[8 pdfad (15+10/111)].^2 13^2 25/4*(17+25/111)^2 26^2]
err_sq_Dan = [0.25*[8 pdfad (15)].^2 13^2 25/4*(17)^2 26^2]
total_Dan = sqrt(sum(err_sq_Dan))
contrib_Dan = 100*err_sq_Dan./sum(err_sq_Dan)

fprintf('\n');
for i=1:length(err_sq_Dan)
   fprintf(1,'%.0f & ', err_sq_Dan(i));   
end
fprintf('\n');

%% Errors with the instrumental uncertainty
fad=0.8; dfad=0.1; pdfad = 100*dfad/fad;
pdfad=30; %From Merk - mean of 0.66 std dev of 0.22 giving approx 30% (actually 33%)
%k uncertainty of 12.5% from Merk too
%17% chosen for re based on VOCALS and LES studies.
%What to choose for tau? Perhaps 17% as for re for now??
re_ins = 10; tau_ins = 10;
%re_ins = 0; tau_ins = 0;%when assuming cancellation

re_err = 17 + re_ins; tau_err = 15 + tau_ins; %sub-pixel plus instrument uncertainty
err_sq_Dan = [0.25*[8 pdfad tau_err].^2 13^2 25/4*re_err^2 30^2]
total_Dan = sqrt(sum(err_sq_Dan))
contrib_Dan = 100*err_sq_Dan./sum(err_sq_Dan)

fprintf('\n');
for i=1:length(err_sq_Dan)
   fprintf(1,'%.0f & ', err_sq_Dan(i));   
end
fprintf('\n');
