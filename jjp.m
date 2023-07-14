%finds closest allowed jjp to target funtion usage: jjp(target)

function jjp(target)
%target=1024;
ms=12; %number of powers to try
f=1.25; %factor of target over which to stop searching

diff=1e39; %big number

for i=0:ms
    for j=0:ms
        for k=0:ms
            prodnew=2^i * 3^j * 5^k;
            diffnew=abs(prodnew-target);
            if diffnew<diff
                diff=diffnew;
                prod=prodnew;
                i1=i;
                i2=j;
                i3=k;
            end
            if (prodnew>target*f) %if gone too far above target
                break
            end
        end
    end
end

fprintf(1,'closest match = %d = 2^%d * 3^%d * 5^%d',prod,i1,i2,i3);
