figure;

for i=25:length(aa)
  r(i).r=( (241-aa(i).i1).^2 + (240-aa(i).i2).^2 ).^0.5;
  vs=repmat(values(i),[size(r(i).r) 1]);
  plot(vs,r(i).r,'kx');
  hold on;
end
