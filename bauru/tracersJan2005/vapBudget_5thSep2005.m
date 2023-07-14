for i=1:length(icediag4)
    iz=2;
    iz2=size(icediag4(i).i,1);
    for j=1:size(icediag4(i).i,3)
        m(i).m(:,j)=TotMassBudgetProf(GridDan(i),icediag4(i).i(:,:,j),GridDan(i).t,iz,iz2);
    end
end