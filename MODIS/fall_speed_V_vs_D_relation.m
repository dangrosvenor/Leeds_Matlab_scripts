function [V]=fall_speed_V_vs_D_relation(D,relation_name)
%LEM / CASIM fall speed relation vs D

if 1==2
%Example usage for plot :-

leg_cell = {'CASIM standard','Abel Shipway 2007 orig uncorrected values','Abel Shipway 2007 just negative a2r','Uplinger'};

D = 10.^ [-5:0.01:-2];
V=fall_speed_V_vs_D_relation(D,leg_cell{1});

figure
plot(D,V); set(gca,'xscale','log'); set(gca,'yscale','log');
set(gca,'ylim',[0.01 20]); set(gca,'xlim',[1e-5 5e-3]);
hold on

V2=fall_speed_V_vs_D_relation(D,leg_cell{2});
plot(D,V2,'b--');

V3=fall_speed_V_vs_D_relation(D,leg_cell{3});
plot(D,V3,'g');

V4=fall_speed_V_vs_D_relation(D,leg_cell{4});
plot(D,V4,'r');

legend(leg_cell);

end



if ~exist('relation_name')
    relation_name = 'Abel Shipway 2007';
end

rho0=1; rho=1; %Ignore this for simplicity

switch relation_name
    case 'CASIM standard'
        %If l_abelshipway is false
        
        a1r = 130;
        b1r = 0.5;
        f1r = 0;
        g1r = 0.5;        
        
        V =   a1r.*D.^b1r .* exp(-f1r.*D) .* (rho0/rho).^g1r;       
        
   case 'Abel Shipway 2007 orig uncorrected values'
         %The wrong ones from the paper
        
        a1r = 4854.1;
        b1r = 1.0;
        f1r = 195;
        g1r = 0.5;                  
        
        a2r = 446;
        b2r = 0.782127;
        f2r = 4085.35;
        g2r = 0.5;
        
        V1 =   a1r.*D.^b1r .* exp(-f1r.*D) .* (rho0/rho).^g1r;
        V2 =   a2r.*D.^b2r .* exp(-f2r.*D) .* (rho0/rho).^g2r;        
        
        V = V1 + V2;   
        
 case 'Abel Shipway 2007 just negative a2r'
         %a2r corrected to be negative (was an erraturm published)
         %Note that the a1r, etc. values are overwritten in the CASIM code
         %in derived_constants.F90 relative to the declarations in
         %mphys_parameteres.F90
        a1r = 4854.1;
        b1r = 1.0;
        f1r = 195;
        g1r = 0.5;                  
        
        a2r = -446;
        b2r = 0.782127;
        f2r = 4085.35;
        g2r = 0.5;
        
        V1 =   a1r.*D.^b1r .* exp(-f1r.*D) .* (rho0/rho).^g1r;
        V2 =   a2r.*D.^b2r .* exp(-f2r.*D) .* (rho0/rho).^g2r;        
        
        V = V1 + V2;        
        
    case 'Uplinger'
        a1r = 4854.1;
        b1r = 1.0;
        f1r = 195;
        g1r = 0.5; 
        
        V =   a1r.*D.^b1r .* exp(-f1r.*D) .* (rho0/rho).^g1r;
        
end

