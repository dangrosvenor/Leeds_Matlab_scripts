function [Nd,Dg,N_aero0,r_single] = activation_run_func(var_in)
%function [Nd,Dg,N_aero,r_single] = activation_run_func(N_aero0,M_aero0,w,r_single,N_type)
% var_in.xxx with xxx=
% N_type - string to describe the type of input
% for Ntype =  'prescribed number'
%  M_aero0  - total aerosol mass MR (kg/kg)
%  N_aero0  - total aerosol number conc (#/kg)
% = 'prescribed single aerosol mass'
%  m_aero_single  - mean aerosol mass (kg)
%  M_aero0
% = 'prescribed median radius (rd)'
%  rd  - median diamter of aerosol dist (m)
%  M_aero0


%Convert all of the variable names in the input structure to actual names
%for ease of use
name_struc='var_in'; %The name of the structure
names = eval(['fieldnames(' name_struc ');']);
for i=1:length(names)
    eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
    eval(eval_str);
end

if ~exist('sig')
    sig = 2; %Lognormal distribution width
end
if ~exist('T')
    T=283; %Temperature in K
end
if ~exist('density')
    density=1777; %aerosol density
end



%% Calculate the aerosol distribution values depending on the type of input
%% values

switch N_type
    case 'prescribed number'
        % N_aero set by function inputs - just calculate r_single for output
        m_aero_single = M_aero0 ./ N_aero0;
        r_single = ( m_aero_single ./ (density .* 4/3*pi) ).^(1/3);
        %Mode diameter in metres (calculated from mass and number)
        Dg = 2* (3.0.*M_aero0.*exp(-4.5.*log(sig).^2.)./(4.0.*N_aero0.*pi.*density)).^(1.0/3.0);
    case 'prescribed single aerosol mass'
%        m_aero_single = density .* 4/3*pi.*r_single.^3;
        r_single = ( m_aero_single ./ (density .* 4/3*pi) ).^(1/3);        
        N_aero0 = M_aero0 ./ m_aero_single;
        %Mode diameter in metres (calculated from mass and number)
        Dg = 2* (3.0.*M_aero0.*exp(-4.5.*log(sig).^2.)./(4.0.*N_aero0.*pi.*density)).^(1.0/3.0);
    case {'prescribed median radius (rd)','prescribed median radius (rd) vary accum'}
        %Supplied M_aero0, but calculate N_aero0
        N_aero0 = (M_aero0.*exp(-4.5.*log(sig).^2.))  ./ (4/3*pi*rd.^3.*density);
        Dg = rd .*2;
        % Not really needed except for outputs:-
        m_aero_single = M_aero0 ./ N_aero0;
        r_single = ( m_aero_single ./ (density .* 4/3*pi) ).^(1/3);
        
end


var_act.N_aero = N_aero0;
var_act.Dg = Dg;
var_act.sig = sig;
var_act.density = density;
var_act.w = w;
var_act.T = T;
var_act.P = P;
var_act.act_routine = act_routine;

switch act_routine
    case {'Nenes','Nenes consistent with Abdul'}
        [Nd] = nenes_func(var_act);  %Function consistent with Adbul is chosen in here
    case 'Abdul'
        [Nd] = abdul_func(var_act);
end





