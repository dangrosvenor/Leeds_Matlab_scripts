
figure('color','w');
leg_str=''

% um_case = 'u-bt398';
% [zNd,Nd_prof] = ukca_bug114_test_v11pt6_FUNC(um_case);
% plot(1e-6*Nd_prof,zNd,'linewidth',8);
% leg_str={'Trunk'};
hold on

um_case = 'u-bt399';
[zNd,Nd_prof,i,j] = ukca_bug114_test_v11pt6_FUNC(um_case);
plot(1e-6*Nd_prof,zNd,'k','linewidth',2);
leg_str=[leg_str {'Test, logical OFF'}];

um_case = 'u-bt400';
[zNd,Nd_prof] = ukca_bug114_test_v11pt6_FUNC(um_case,0,i,j);
plot(1e-6*Nd_prof,zNd,'b--','linewidth',2);
leg_str=[leg_str {'Test, logical ON'}];

legend(leg_str);
ylabel('Height (m)');
xlabel('CDNC (cm^{-3})');

savename = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case '/ukca_bug114_plots'];
%saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);

