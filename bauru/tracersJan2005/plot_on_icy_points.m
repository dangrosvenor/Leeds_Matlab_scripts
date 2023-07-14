figure
pcolor(Grid.Y1,Grid.Z+620,sum(TwoD.Q(:,:,[10]),3));
pcolor(Grid.Y1,Grid.Z+620,sum(TwoD.Q(:,:,[4:6]),3));

shading interp;
colorbar;

hold on;
for k=1:250
	ii=find(sum(TwoD.Q(k,:,[7]),3)>1e6);
    ii=find(sum(TwoD.Q(k,:,[10]),3)>-99 & TwoD.W(k,:)>0.5 & sum(TwoD.Q(k,:,[2 6]),3)>1e-7);
	ii=find(sum(TwoD.Q(k,:,[4:6]),3)>5e-3);

	for i=1:length(ii)
		plot(Grid.Y1(ii(i)),Grid.Z(k)+620,'wx');
	end
end