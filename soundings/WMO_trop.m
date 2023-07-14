function [trop,trop_press,trop2,trop_press2]=wmo_trop(H,temp2,press,smooth,windowSize)
%function [trop,trop_press,trop2,trop_press2]=wmo_trop(H,temp2,press,smooth,windowSize)
%calculates the WMO definition of the tropopause
%using fields of temp, heights and press


%definition states:

%The first tropopause (i.e., the conventional tropopause) is defined as the 
%lowest level at which the lapse rate decreases to 2 K/km or less, and the average lapse rate from this 
%level to any level within the next higher 2 km does not exceed 2 K/km. 

%If above the first tropopause the average lapse rate between any level and 
%all higher levels within 1 km exceed 3 K/km, then a second tropopause is defined 
%by the same criterion as under the statement above. This tropopause may be either within or above the 1 km layer. 

%A level otherwise satisfying the definition of tropopause, but occuring at 
%an altitude below that of the 500 mb level will not be designated a tropopause unless it is 
%the only level satisfying the definition and the average lapse rate fails to exceed 3 K/km 
%over at least 1 km in any higher layer.
%

if smooth==1
%smooth the temperature data for easier viewing of lapse rates
	temp2=filter(ones(1,windowSize)/windowSize,1,temp2);
end


%calculate lapse rate in K/km
LR=1000*diff(temp2)./diff(H);

i2=find(LR>-2);  %find points where exceeds 2 K/km

for i=1:length(i2)
    trop_found=1;
    min_mean=9e99;
    if (H(end)-H(i2(i))>2100) & (press(i2(i))<500e2) %make sure is below 500mb and have enough points for 2 km
        i2000=findheight(H,H(i2(i))+2000); %find index 2 km above
        for h=i2(i):i2000 %doing running averages over all points from height concerned up to 2 km
            mean_lap=mean(LR(i2(i):h));
            min_mean=min(min_mean,mean_lap);
            if mean_lap<-2   %if mean LR is not above 2 K/km then not tropopause
                trop_found=0;
                break
            end
            
        end
        if trop_found==1
            itrop=i2(i);            
            trop=H(itrop); trop_press=press(itrop)/100;
            break
        end
    end
end

%look for second trop
found=0;
for i=itrop:length(H)
    trop_found=1;
    min_mean2=9e99;
    if (H(end)-H(i)>1100) & (press(i)<500e2)
        i1000=findheight(H,H(i)+1000);
        for h=i:i1000
            mean_lap=mean(LR(i:h));
            min_mean=min(min_mean,mean_lap);
            if mean_lap>-3   %average up to ALL levels within  1km must exceed 3 K/km
                trop_found=0;
                break
            else
                '';
            end
            
        end
        if trop_found==1
            trop2=H(i); trop_press2=press(i)/100;
            found=1;
            break
        end
    end
end
if found==0; 
    trop2=0;
    trop_press2=0;
    fprintf(1,'\nNo second tropopause found.'); 
end