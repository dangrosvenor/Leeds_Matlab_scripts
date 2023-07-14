dth_thresh = 0.015;
hrange = 1000; %distance (m) over which condition must be satisfied
hmin=7e3; %minimum height to look at
% 
% for i=1:size(TwoDf.TH2,2)
%     th=TwoD.TH2(:,i);
%     dth=diff(th)./diff(GridDan(1).Z);
%     ii=find(dth>dth_thresh);
%     iibl=find(GridDan(idir).Z(ii)+620 > 7e3);
%     htrop(i)=GridDan(1).Z( min(ii(iibl)) ) + 620;
% end
% 
% 'done'

for i=1:size(TwoD.TH2,2)
    th=TwoD.TH2(:,i);
    [t,htrop(i)]=tropopause_th(th,GridDan(idir).Z+620,dth_thresh,hrange,hmin);
end

'done'
