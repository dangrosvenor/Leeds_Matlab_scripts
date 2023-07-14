dirpath='C:\Documents and Settings\Login\My Documents\samis_model\icemodel1\';

% out_folder='output_1st_run\';
% out_folder='output\';
% out_folder='output_HM_20\';


figure;
bins=30;
bin=60;
bin2=60;

iplot=3;


for III=1:1:1
  hold off
  if(III<10)
    figure;
    
%    eval(['load ' dirpath out_folder 'output00' num2str(III) '.out'])
%    eval(['load ' dirpath out_folder num2str(III) '.out'])


% 
%     fname_out=['output00' num2str(III)];
%     load([dirpath out_folder fname_out '.out'])
%     fname_ice=['output00' num2str(III)];
%     load([dirpath out_folder fname_ice '.out'])

%     a=eval(fname_ice);
%     b=eval(fname_out);

       % cmodel_read_and_plot % - uses this now to plot - this loads in a and the variables seperately
        
    a=aa(iplot).dat;    
     
    ax1=meshc(log10(a(1,2:bin+1)),a(:,1),log10((a(:,bin+3:2*bin+2)+1e-6)./ ...
                                        (log10(a(:,3:bin+2))-log10(a(:,2:bin+1)))));
    axis([-7 -2.5 0 4700 -7 1])
    view(2)
    set(ax1(1),'FaceColor','flat')
    [sum(a(1:72,246:249),1) sum(sum(a(1:72,246:249),1))] 
  elseif(III<100)
     figure(10)
    
     eval(['load output/iceoutf0' num2str(III) '.out'])
     eval(['load output/outputf0' num2str(III) '.out'])
     eval(['a=iceoutf0' num2str(III) ';'])
     eval(['b=outputf0' num2str(III) ';'])
     ax1=meshc(log10(a(1,2:bin+1)),a(:,1),log10((a(:,bin+3:2*bin+2)+1e-6)./ ...
                                               (log10(a(:,3:bin+2))-log10(a(:,2:bin+1)))));
     
    axis([-7 -2.5 0 4700 -7 1])
    view(2)
    set(ax1(1),'FaceColor','flat')
    [sum(a(1:72,246:249),1) sum(sum(a(1:72,246:249),1))] 
  else
    figure(10)
    
    eval(['load output/iceoutf' num2str(III) '.out'])
    eval(['load output/outputf' num2str(III) '.out'])
    eval(['a=iceoutf' num2str(III) ';'])
      eval(['b=outputf' num2str(III) ';'])
     ax1=meshc(log10(a(1,2:bin+1)),a(:,1),log10((a(:,bin+3:2*bin+2)+1e-6)./ ...
                                               (log10(a(:,3:bin+2))-log10(a(:,2:bin+1)))));
     
     axis([-7 -2.5 0 4700 -4 1])
     view(2)
    set(ax1(1),'FaceColor','flat')
    [sum(a(1:72,246:249),1) sum(sum(a(1:72,246:249),1))] 
  end 
  III
  pause(.5)
  
  figure;
  ax1=meshc(log10(a(1,2*bin+3:2*bin+2+bin2)),a(:,1),log10((a(:,2*bin+4+bin2 ...
                 :2*bin+2*bin2+3)+1e-6)./(log10(a(:,2*bin+4:2*bin+bin2+3))-log10(a(:,2*bin+3:2*bin+bin2+2)))));
  axis([-8 -3 0 4700 -7 1])
     
  view(2)
  set(ax1(1),'FaceColor','flat')
  hold on
  semilogx(log10(b(:,2:bins+1)),b(:,1),'m')
end

figure
plot(time,c_ice_tot);