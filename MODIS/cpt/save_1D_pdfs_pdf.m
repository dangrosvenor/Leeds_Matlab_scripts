
switch x_axis_vals
    case {'Re COSP GCM','CDR Polder2'}

        var_str = 'Re';
        
end


eval(['Xbins_' var_str '_' gcm_str '= xdat(1).x;']);
eval(['PDF_' var_str '_' gcm_str '= ydat(1).y;']);

['save(' var_str '1Dpdfs2.mat,PDF_' var_str '_' gcm_str ',Xbins_' var_str '_' gcm_str ',''-APPEND'')'];
