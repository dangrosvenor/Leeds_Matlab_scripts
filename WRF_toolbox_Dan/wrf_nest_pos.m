function ix=wrf_nest_pos(dom,dom2,ratio)
% function ix=wrf_nest_pos(dom,dom2,ratio)
% works out the position of domain 2 (bottom left) for given number of 
% points for domain one and two and grid ratio
% assuming that the nested domain is in the middle.
ix = dom/2 - dom2/2/ratio;