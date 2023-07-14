function [mids]=mid_vals(X)
%function [mids]=mid_vals(X)
%Finds the midpoints of a 1D array - mids=0.5*(X(1:end-1)+X(2:end));
mids=0.5*(X(1:end-1)+X(2:end));