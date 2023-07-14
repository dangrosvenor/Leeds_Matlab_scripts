function hr=heart_rate(e,restHR,maxHR)

% e=65; %percentage effort
% maxHR=205;
% restHR=55;

hr=e/100*(maxHR-restHR) + restHR;