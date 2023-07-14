function tdew=Tdew(Q,P)        
%function tdew=Tdew(Q,P)  
%Q is the vapour MR in kg/kg
%P is the air pressure in Pa
%Uses the Goff-Gratch formula (GGEW function) - note, this is the
%same as one used in SatVapPress with 'goff' selcted

%if(length(Q)>1)
 %   fprintf('Error\n');
    %else
    
    T=273.15; %this is just a starting point for the fzero iteration
    f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
        

    
  
                for i=1:length(Q)
                    tdew(i)=fzero(inline(sprintf('q_sat(x,%f)-%.20f;',P(i),Q(i))),T);
                end
                %so this runs solves the function y=q_sat(x,P(i))-Q(i) where x is the temperature that fzero is trying to solve for
                %and finds y=0, i.e. when Qsat=Q(i)

        

    
    %end
