function qv = vap_from_Tdew(Tdew,P)
%function qv = vap_from_Tdew(Tdew,P)
%returns the vapour mixing ratio in kg/kg for a given dew point temperature (K) and air pressure (Pa)

f=1e6*28.97/18;

qv = satvappress(Tdew,'goff','liq',P,1);

qv = qv/f;