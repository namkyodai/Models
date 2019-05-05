SV=dlmread('SV.txt');
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
   for f=0:c-1
        sv=sum(sv0(1:i))+f*(sum(sv1(1:i))-sum(sv0(1:i)))/c;
        if sv<P&j==0
            le(1+i)=BB-1+a*(f)/c;
            d=1;
            j=1;
        end
%         if f/c==1
%             d=1;
%         end
%         f=f+1;
    end
end
end
le
