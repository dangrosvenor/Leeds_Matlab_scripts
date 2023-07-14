function [out] = load_mat_to_var(filename,var)

out = load(filename,var);
eval(['out = out.' var ';']);