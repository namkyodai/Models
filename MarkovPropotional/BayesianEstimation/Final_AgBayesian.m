clear all
% HHHH=-1;
% while HHHH<0
A=1;%time interval
%% Sample generation
%%% Initial setting
%state classification
I=5;
Si=[0,0.24,0.48,0.48,0.48];
Sp=[0,0,0,0.25,0.5];
%%% Initiation
%const.
d=20;
Cs=2.4;
%parameters
meanDcl=10;
varDcl=10;
% Ccl=Cs*(1-erf(d/(2*sqrt(Dcl*t))));


%%% Propagation
%const.
%parameters
meanXpro=0.0004;
varXpro=0.00001;
meanw0=0.05;
varw0=0.001;

% t=10;
% w0=0.05;Xpro=0.004;
% w=w0+t^2*Xpro/2;

% simulation
N=10000;
Z=100;
% A=1;% time interval
iDB=zeros(Z/A,N);
pDB=zeros(Z/A,N);
for n=1:N
     Dcl=-1;
        while Dcl<0
        Dcl=randn*sqrt(varDcl)+meanDcl;
        end
        w0=-1;
        while w0<0
            w0=randn*sqrt(varw0)+meanw0;
        end
        Xpro=-1;
        while Xpro<0
            Xpro=randn*sqrt(varXpro)+meanXpro;
        end
%     for z=A:A:Z
%         JJJ=-1;
%         while JJJ<0
%         t=A*(z-A);
%        
% 
%         iDB(round(z/A),n)=Cs*(1-erf(d/(2*sqrt(Dcl*t))));
%       
%         pDB(round(z/A),n)=w0+t^2*Xpro/2;
%         if z==A
%             if iDB(round(z/A),n)>=0 && pDB(round(z/A),n)>=0
%                 JJJ=1;
%             end
%         else
%             if iDB(round(z/A),n)>=0 && pDB(round(z/A),n)>=0 && iDB(round(z/A),n)>=iDB(round(z/A)-1,n) && pDB(round(z/A),n)>=pDB(round(z/A)-1,n)
%                 JJJ=1;
%             end
%         end
%         end
%     end
    for z=A:A:Z
%         JJJ=-1;
%         while JJJ<0
        t=z;
       

        iDB(round(z/A),n)=Cs*(1-erf(d/(2*sqrt(Dcl*t))));
      
        pDB(round(z/A),n)=w0+t^2*Xpro/2;
%         if z==A
%             if iDB(round(z/A),n)>=0 && pDB(round(z/A),n)>=0
%                 JJJ=1;
%             end
%         else
%             if iDB(round(z/A),n)>=0 && pDB(round(z/A),n)>=0 && iDB(round(z/A),n)>=iDB(round(z/A)-1,n) && pDB(round(z/A),n)>=pDB(round(z/A)-1,n)
%                 JJJ=1;
%             end
%         end
%         end
    end
end




%% Bayesian estimation

DB=zeros(Z/A,N);
DB(1,:)=ones(1,N);
for n=1:N
    for z=1:Z/A-1
%     for ji=1:I-1

    if 0<=iDB(z,n) && 0.24>=iDB(z,n)
    DB(z+1,n)=1;
    elseif 0.24<iDB(z,n)  && 0.48>=iDB(z,n)
        DB(z+1,n)=2;
    elseif 0.48<iDB(z,n) && 0<=pDB(z,n) && 0.25>=pDB(z,n)
        DB(z+1,n)=3;        
    elseif 0.48<iDB(z,n) && 0.25<pDB(z,n) && 0.5>=pDB(z,n)
        DB(z+1,n)=4; 
    end

%     end
    if DB(z,n)==0;
DB(z,n)=I;
    end
    end
  if DB(end,n)==0;
DB(end,n)=I;
    end
end
        Prop=zeros(Z/A,I);
for z=1:Z/A
Prop(z,:)=[length(find(DB(z,:)==1)),length(find(DB(z,:)==2)),length(find(DB(z,:)==3)),length(find(DB(z,:)==4)),length(find(DB(z,:)==5))]./N;
end
            

%% Estimation of Markov transition probability


% A=round(A/3)
% Initial setting

alpha0=0.001.*[0.8,0.0529,0.0529,0.0529,0.0529;
    0,0.96,0.03,0.01,0.01;
    0,0,0.93,0.1,0.1;
    0,0,0,0.93,0.1];% Parameter of prior distribution (Dirichlet distribution)

% alpha0(1,:)=alpha0(1,:)./2;
%alpha02=[5,5,5,5];
%alpha03=[5,5,5];
%alpha04=[5,5];
G1=1.5;G2=3;% hyper parameter
% P0=[0.5,0.1,0.1,0.05,0.05;
%     0,0.5,0.2,0.2,0.1;
%     0,0,0.7,0.2,0.1;
%     0,0,0,0.9,0.1;
%     0,0,0,0,1];% initial MTP
P0=[    0.7617    0.2278    0.0061    0.0022    0.0022
         0    0.8354    0.1448    0.0143    0.0055
         0         0    0.8691    0.1018    0.0291
         0         0         0    0.8581    0.1419
         0         0         0         0    1];% initial MTP

Sam=11000;% # of roops of MCMC
samP=zeros(Sam,I*I);% samples is recorded in this matrix
samalpha=zeros(Sam,I*(I-1));
Rej=zeros(1,I-1);% initial rejection rate

% MCMC
P=P0;%alpha=alpha0;
for sam=1:Sam
     tic
    for tP=1:I-1
        for ttP=tP:I-1
      
%         for ta=1:I-1
    L=0;% log likelihood
    for a=1:Z/A-1
        for n=1:N
        l=P(DB(a,n),DB(a+1,n));
%         end
                
            
        
            
%         l=l/(1-kl);    
        L=L+log(l);% Natural logarithm 
        end
    end
    
    %Prior distribution
    Prior=1;
    for ddd=tP:I
    Prior=Prior*(P(tP,ddd)^(alpha0(tP,ddd)-1))/gamma(alpha0(tP,ddd));
    end
    Prior=Prior*gamma(sum(alpha0(tP,tP:I)));
       L=L+log(Prior);
    

    NP=P;
%     for i=tP:tP
%         Y=10000;
%         RG=zeros(Y,I-i+1);
%         for u=1:I-i+1
%             RG(:,u)=gamrnd(alpha(i,i+u-1),1,Y,1);
%         end
%         for u=1:I-i+1
%         NP(i,i+u-1)=sum(RG(:,u))/sum(sum(RG));
%         end
%     end
%     NP(end,end)=1;
JJJ=2;
while JJJ>1
NP(tP,ttP)=P(tP,ttP)+randn/500;%10;
if sum(NP(tP,:))-NP(tP,ttP+1)>1||sum(NP(tP,:))-NP(tP,ttP+1)<0
    JJJ=2;
elseif NP(tP,ttP)>1||NP(tP,ttP)<0
    JJJ=2;
else
    JJJ=0;
end
end

NP(tP,ttP+1)=1-(sum(NP(tP,:))-NP(tP,ttP+1));
   
        NL=0;% log likelihood
   for a=1:Z/A-1
        for n=1:N
        Nl=NP(DB(a,n),DB(a+1,n));
%         end
                
            
        
            
%         l=l/(1-kl);    
        NL=NL+log(Nl);% Natural logarithm 
        end
   end
    
       %Prior distribution
    NPrior=1;
    for ddd=tP:I
    NPrior=NPrior*(NP(tP,ddd)^(alpha0(tP,ddd)-1))/gamma(alpha0(tP,ddd));
    end
    NPrior=NPrior*gamma(sum(alpha0(tP,tP:I)));
       NL=NL+log(NPrior);
%         for ui=1:I-1
%             S=1;
%             for ab=ui:I
%                 S=S*gamma(Nalpha(ui,ab));
%             end
%             D=S/gamma(sum(Nalpha(ui,:)));
%             di=1;
%             for uii=ui:I-1
%                 di=di*NP(ui,uii)^(Nalpha(ui,uii));
%             end
%         NL=NL+log(di/D);
%         end
%    for ai=1:I-1
%         for aii=1:I
%             if Nalpha(ai,aii)==0
%             else
%             NL=NL+log(gampdf(Nalpha(ai,aii),G1,G2));
%             end
%         end
%    end
jpre=normpdf((1-P(tP,ttP))/(sqrt(500)))/(sqrt(500))/(normcdf((1-P(tP,ttP))/(sqrt(500)))-normcdf((0-P(tP,ttP))/(sqrt(500))));
jcan=normpdf((1-NP(tP,ttP))/(sqrt(500)))/(sqrt(500))/(normcdf((1-NP(tP,ttP))/(sqrt(500)))-normcdf((0-NP(tP,ttP))/(sqrt(500))));

mjudge=min(1,exp(NL+log(jpre)-L-log(jcan)));

        ur=rand();
        if ur > mjudge
             Rej(tP)=Rej(tP)+1;
        else
           P=NP;%alpha=Nalpha;
        end
%         end
% toc
% sam

% tP
% P
        end
    end
    
%  
samP(sam,:)=reshape(P',[1,I*I]);
    
        sam
        
        toc
        P
        Rej./sam/(I-1)
end

%% Credit interval
LE=zeros(1,Sam);
for samm=1:Sam
    St=zeros(1,I);St(1)=1;
    Pi=reshape(samP(samm,:),[I,I])';
    for tz=1:A:1000000
        St=St*Pi;
        if St(end)>0.5
            LE(samm)=tz-A;
            break
        end
    end
end
    

%     0.8882    0.1114    0.0001    0.0001    0.0001
%          0    0.8188    0.0953    0.0512    0.0347
%          0         0    0.9155    0.0844    0.0002
%          0         0         0    0.8810    0.1190
%          0         0         0         0    1.0000
         
%              0.7888    0.1902    0.0109    0.0060    0.0041
%          0    0.6704    0.1653    0.0950    0.0692
%          0         0    0.8381    0.1516    0.0104
%          0         0         0    0.7762    0.2238
%          0         0         0         0    1.0000

%    0.7926    0.2073    0.0000    0.0000    0.0000         0
% 
%   7 —ñ‚©‚ç 12 —ñ
% 
%     0.6370    0.1758    0.1111    0.0761         0         0
% 
%   13 —ñ‚©‚ç 18 —ñ
% 
%     0.8476    0.1524    0.0000         0         0         0
% 
%   19 —ñ‚©‚ç 24 —ñ
% 
%     0.7674    0.2326         0         0         0         0
% 
%   25 —ñ
         