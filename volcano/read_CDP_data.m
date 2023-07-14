time_CDP=nc{'Time'}(:);

CDP_bins = [2.0000000, 3.0000000, 4.0000000, 5.0000000, 6.0000000,...
                        7.0000000, 8.0000000, 9.0000000, 10.000000,...
                        11.000000, 12.000000, 13.000000, 14.000000, 16.000000, 18.000000,...
                        20.000000, 22.000000, 24.000000, ...
                        26.000000, 28.000000, 30.000000, 32.000000, 34.000000, 36.000000, ...
                        38.000000, 40.000000, 42.000000, ...
                        44.000000, 46.000000, 48.000000, 50.000000];
                    
                    

clear data_particle
                    data_particle(:,1)=zeros(size(time_CDP))';
                    
                    for ichannel=1:30    
                        if ichannel<10
                            channel_str=['0' num2str(ichannel)];
                        else
                            channel_str=[num2str(ichannel)];
                        end
                        
                        eval_str = ['data_particle(:,ichannel+1) = nc{''CDP_' channel_str '''}(:);'];
                        eval(eval_str);
                    end
                    
                    