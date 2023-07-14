function T=tempLES(Grid)
%function T=tempLES(Grid)

T=Grid.THREF./(1e5./Grid.PREFN).^0.286;