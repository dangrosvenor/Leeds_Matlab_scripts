itext_in_box = 1; 

% box_region = '1'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
% %box_region = '2'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
% box_region = '3'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
box_region = '4'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map %N. Atlantic basin
% %box_region = '5'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
% %box_region = '6'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
% %box_region = '7'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
% box_region = '8'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '9'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '10'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '11'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '12'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '13'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '14'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '15'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '17'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map

%box_region = '21'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '22'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '23'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map

%box_region = '06'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '07'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '08'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '09'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = '10'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = 'L01'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = 'L02'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = 'L03'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = 'L04'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map
%box_region = 'L05'; ACSIS_Robson_paper_choose_regional_box2; ioverride_box_colour=1;plot_box_on_map

%Global paper regions
% Oceans - 00 (global 60S to 60N), Antarctic sea-ice region (04), Arctic sea-ice region (05), 
%  SO (03), Southern Pacific (06),  South Indian Ocean (07), S. Atlantic, Tropical Pacific, Tropical Atlantic, Tropical Indian Ocean, Northern Pacific, Northern Atlantic
% Land - South America (L01), N. America (L02), Africa (L03),  Europe+N.
% Asia (L04), S. Asia (L05).

%ACtually, simplfy to -40 to 0 and 0 to +40 for oceans.
%Global paper regions
% Oceans - 00 (global 60S to 60N), Antarctic sea-ice region (04), Arctic sea-ice region (05), 
%  SO (03), Southern Pacific (06),  Indian Ocean (07), S. Atlantic (08),
%  Northern Atlantic (09), Northern Pacific (10)
% Land - South America (L01), N. America (L02), Africa, Europe+Asia.


%reset back to original value
if exist('box_region_DRIVER')
    box_region = box_region_DRIVER; ACSIS_Robson_paper_choose_regional_box2;
end
