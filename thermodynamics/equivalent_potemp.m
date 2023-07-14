function equiv=equivalent_potemp(TK,P,qv)
%function equiv=equivalent_potemp(TK,P,qv)
f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

L=size(TK,1);
for i=1:size(TK,1)
    
    if size(P,1)~=size(TK,1)
        fprintf(1,'\n');
        disp('*****  Problem with size of P ********');
        fprintf(1,'\n');
        return
    end

    
    if size(qv,1)~=size(TK,1)
        fprintf(1,'\n');
        disp('*****  Problem with size of qv ********');
        fprintf(1,'\n');
        return
    end

%AMS definition - is pretty similar
qsat=SatVapPress(TK(i,:),'goff','liq',P(i,:),1)/f;   
RH = qv(i,:)./qsat;
equiv(i,:) = TK(i,:).*(1e5./P(i,:)).^0.286 .* RH.^(-qv(i,:)*4.615e2/1005.7)...
    .* exp( 2.501e6*qv(i,:)./TK(i,:)/1005.7 );

end
                        