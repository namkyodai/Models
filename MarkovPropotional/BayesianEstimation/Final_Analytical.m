clear all
a=1;%time interval
%% Sample generation
%%% Initial setting
%state classification
% I=5;
Si=[0,0.24,0.48,0.48,0.48];
%%% Initiation
%const.
d=20;
Cs=2.4;
%parameters
meanx=10;
varx=10;
minx=0;
tZ=100;

% simulation
I=1;J=1;
N=1000000;
P=0;p1=0;p2=0;
    V=rand(N,1);F1=0;F2=0;F=0;
for n=1:N
    v0=V(n);
      tz=tZ*v0;
            overx=1/(tz+a)*(d/(2*erfinv(1-0.24/(Cs))))^2;
        overy=min(0.24,Cs*(1-erf(d/(2*sqrt(overx*tz)))));
        x1=1/(tz)*(d/(2*erfinv(1-overy/(Cs))))^2;
        delta1=1-normcdf((minx-meanx)/sqrt(varx));
        f1=(normcdf((x1-meanx)/sqrt(varx))-normcdf((minx-meanx)/sqrt(varx)))/delta1;
        F1=F1+f1/tZ;
       
        x2=1/(tz)*(d/(2*erfinv(1-0.24/(Cs))))^2;
        delta2=1-normcdf((minx-meanx)/sqrt(varx));
        f2=(normcdf((x2-meanx)/sqrt(varx))-normcdf((minx-meanx)/sqrt(varx)))/delta2;
        F2=F2+f2/tZ;
        F=F+f1/f2;
        f1/f2;
        tz;
 if mod(n,1000)==0
    n
 end
end
FF1=F1/N*tZ
FF2=F2/N*tZ
F1/F2;
PI=FF1/FF2



