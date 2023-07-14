function w=ws(t,p)

epsilon=0.622;
es=SatVapPress(t,'lem','ice',p);
w=epsilon*es./(p-es);