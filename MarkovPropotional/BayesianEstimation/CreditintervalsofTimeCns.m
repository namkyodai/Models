AcLE=zeros(1,5);

for eee=1:1
      SV=zeros(500,5);ssv=[1 0 0 0 0]; %Pi=reshape(samP(eee,:),[I,I])';
      Pi=[0.0296	0.0973	0.2036	0.2273	0.4421
0	0.0227	0.0619	0.0666	0.8488
0	0	0.323	0.3655	0.3114
0	0	0	0.0259	0.9741
0	0	0	0	1];
    for tt=1:500
        ssv=ssv*Pi;SV(tt,:)=ssv;
    end

[bb,b]=size(SV);
le=zeros(1,5)';
a=1;
c=1000;
P=0.5;
for i=1:4 
    d=0;j=0;
for BB=1:bb
    
    if BB==1
        sv0=[1,0,0,0,0];
        sv1=SV(1,:);
    else
        sv0=SV(BB-1,:);sv1=SV(BB,:);
    end
%    for f=0:c-1
        sv=sum(sv0(1:i));%+f*(sum(sv1(1:i))-sum(sv0(1:i)))/c;
        sssv=sum(sv1(1:i));
        if sv>P && sssv<=P && j==0
            le(1+i)=BB-1+a*(sv-P)/(sv-sssv);
            d=1;
            j=1;
            
        end
%         if f/c==1
%             d=1;
%         end
%         f=f+1;
%     end
end
end
AcLE(eee,:)=le';
eee
end
AcLE
% CreLE=zeros(1,2);
% AccLE=sort(AcLE(:,end));
% CreLE=[AccLE(10000*0.05),AccLE(10000*0.95)]
