%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%     program name : Markov.m  (matlab program)                     %
%                                                                   %
%                                 2009.10  Ryosuke OKIZUKA          %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   このプログラムは，マルコフモデル(exp付き)のパラメータを推計します．
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%% 初期設定 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileinp1=['sample.txt'];      %-- 入力ファイル名

In=6;                          %-- 入力ファイル列数

Maxi=10;                      %-- ニュートンラフソン法の最大繰返数
Eps=1e-06;                     %-- ニュートンラフソン法の収束条件

%%--1--before
%%--2--after
%%--3--interval


Cd=[5 6];                      %-- 説明変数として採用する変数　上参照．
XBP=[1 1 1];                   %-- 採用変数
ZP=2;                          %-- 採用経過年　　　　　　　
Jmax=7;                        %--段階数



XB=[-1.0401 0 3.15 -1.486 0 3.33 -1.957 0.7166 0 -2.4399 0.8705 0.5148 -2.3599 0 0 -1.9984 1.5473 0]';   %--初期パラメータ
XP=[1 3 4 6 7 8 10 11 12 13 16 17];    %--推計するパラメータの列

DE=[2 1; 2 2; 3 3; 2 5; 3 5; 3 6];            %--削除パラメータ(変数，段階)
DL=[2 5 9 14 15 18];                       %--削除パラメータ（XBにおける列数）

DDLL=0;                                  %--削除しているかどうか？ 1--削除無し

xp=length(XP);
[md ml]=size(DE);
pl=length(XB);

%%%%%%%% ファイルの読み込み %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp1=fopen(fileinp1,'r');
DB=fscanf(inp1,'%f',[In,inf]);

%%%%%%%% パラメータ設定 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Be=DB(1,:);                  % Before
Af=DB(2,:);                  % After
Ins=DB(3,:);                 % Interval


pk=length(Cd);
for i=1:pk
     X(i,:)=DB(Cd(i),:);     % Explanation Value
end

[M,K]=size(X);               %-- サンプルサイズ：Kの数え上げ

X1=ones(1,K);                %-- 定数項の設定
rXk=[X1;X];                  %-- 定数項＋その他の変数  

%%%%%%% 最大値１となるように規格化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:M+1
         Xk(i,:)=rXk(i,:)/max(rXk(i,:));
    end

%%%%%%% 大ループの開始 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for ii=1:Maxi
    GR=zeros(pl,1);
    HE=zeros(pl,pl);
    
%%%%%%% 中ループの開始 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for k=1:K
     
%%%%%%% 状態の把握 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=Be(k);
j=Af(k);

%%%%%%% θi(i=l,n)の値 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   BPX=zeros(1,Jmax);
   Theta=zeros(1,Jmax);            
             
    for l=1:Jmax-1         
        for h=1:pk+1 
            BPX(l)=BPX(l)+XB((l-1)*(pk+1)+h)*Xk(h,k);
        end
    end
    Theta(1:Jmax-1)=double(exp(BPX(1:Jmax-1)));
           
%%%%%%% 各θi-θjの計算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
%%--対各要素,サンプルデータの範囲外は，1に置換

Thetasa=ones(Jmax,Jmax);

  for s=i:j
      for q=i:j
        if s~=q
        Thetasa(s,q)=Theta(s)-Theta(q);    
        end
      end
  end

%%%%%%% 一階偏微分の計算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df=zeros(1,Jmax-1);

  for l=i:j
    if l~=Jmax
       sum2=0;                   %%--分母
       for h=i:j
           sum1=0;               %%--分子
           if h==l
              for p=i:j
                  if p~=l
                     sum1=sum1+1/(Thetasa(p,l));
                  end
              end
              sum1=sum1-Ins(k);
           else
           sum1=sum1-1/Thetasa(l,h);
           end
           if l~=j
              sum1=sum1+1/Theta(l);
           end
           RESERVE=double(exp(-Theta(h)*Ins(k)))/prod(Thetasa(:,h));
           
           sum1=sum1*RESERVE;
           sum2=sum2+RESERVE;
           df(l)=df(l)+sum1;
       end
       df(l)=df(l)/sum2;
    end
  end
                      
%%%%%%% 二階偏微分の計算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HES=zeros(Jmax-1,Jmax-1);

%%--非対角要素

  for l=i:j-1
    for n=l+1:j
        if n~=Jmax
           sum2=0;               %%--分母
           for h=i:j
               sum1=0;           %%--分子
               if h==l
                  for p=i:j
                      if p~=l
                         sum1=sum1+1/Thetasa(p,l);
                      end
                  end
                  sum1=sum1+1/Theta(l)+1/Thetasa(n,l)-Ins(k);
                  sum1=-sum1/Thetasa(n,l);
               elseif h==n
                      for p=i:j
                          if p~=n
                              sum1=sum1+1/Thetasa(p,n);
                          end
                      end
                      sum1=(sum1-Ins(k))*((1/Theta(l))-(1/Thetasa(l,n)));
                      sum1=sum1-(1/((Thetasa(l,n))^2));
               else 
                      sum1=(1/Theta(l))-(1/Thetasa(l,h));
                      sum1=-sum1/Thetasa(n,h);
               end
               RESERVE=(exp(-Theta(h)*Ins(k)))/prod(Thetasa(:,h));
               HES(l,n)=HES(l,n)+sum1*RESERVE;
               sum2=sum2+RESERVE;
           end
           HES(l,n)=-df(l)*df(n)+HES(l,n)/sum2;
           if n~=j
              HES(l,n)=HES(l,n)+df(l)/Theta(n);
           end
           HES(n,l)=HES(l,n);
        end
    end
  end

%%--対角要素

  for l=i:j
    if l~=Jmax
       sum2=0;                 %%--分母
       for h=i:j
           sum1=0;             %%--分子
           for p=i:j
               if p~=l
                  sum1=sum1+1/((Thetasa(p,l))^2);
               end
           end
           RESERVE=-Ins(k);
           if h==l
              for p=i:j
                  if p~=l
                     RESERVE=RESERVE+1/Thetasa(p,l);
                  end
              end
              sum1=sum1+(RESERVE)^2;
              if l~=j
                 sum1=sum1+(RESERVE-(1/Theta(l)))/Theta(l);
              end
           else
              sum1=2/((Thetasa(l,h) )^2);
              if l~=j
                 sum1=sum1-((1/Theta(l))^2);
                 sum1=sum1-1/(Theta(l)*Thetasa(l,h));
              end
           end
           RESERVE=double(exp(-Theta(h)*Ins(k)))/prod(Thetasa(:,h));
               HES(l,l)=HES(l,l)+sum1*RESERVE;
               sum2=sum2+RESERVE;
       end
       HES(l,l)=-((df(l))^2)+HES(l,l)/sum2;
           if l~=j
              HES(l,l)=HES(l,l)+df(l)/Theta(l);
           end
    end
  end
  
%%%%%%% 削除パラメータの算定 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kX=zeros(pk+1,Jmax-1);

  for l=i:j
    for M=1:pk+1
        kX(M,l)=Xk(M,k);
    end
  end
  for dl=1:md
      kX(DE(dl,1),DE(dl,2))=0;
  end
%%%%%%% グラディアントベクトルの計算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   for l=i:j
       if l~=Jmax
          for M=1:pk+1
              GR((l-1)*(pk+1)+M)=GR((l-1)*(pk+1)+M)+df(l)*kX(M,l)*Theta(l);
          end
       end
   end
       
%%%%%%% ヘシアン行列の計算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for l=i:j
      for n=i:j
          for M=1:pk+1
              for M2=1:pk+1
                  if l~=Jmax
                     if n~=Jmax
                        HE((l-1)*(pk+1)+M,(n-1)*(pk+1)+M2)=HE((l-1)*(pk+1)+M,(n-1)*(pk+1)+M2)...
                          +HES(l,n)*kX(M,l)*kX(M2,n)*Theta(l)*Theta(n);
                     end
                  end
              end
          end
      end
  end
 end                            %%--中ループ終了

%%%%%%% 削除パラメータ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if DDLL~=1
      GR(DL)=[];
      HE(:,DL)=[];
      HE(DL,:)=[];
      XB(DL)=[];
   end
   
%%%%%%% 更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 MXB=XB-inv(HE)*GR;                         %-- 非線形連立方程式の直接解法

 Er(ii)=max(double(abs((MXB-XB)./MXB)));
                                            %-- 収束判定指標
    if Er(ii) < Eps                         %-- 収束状況の確認
         break
    else
       xx=0;        
       XB=zeros(pl,1); 
       for xd=1:xp
           xx=xx+1;
           XB(XP(xd))=MXB(xx);
       end      
    end
end
  
%%%%%%% 推計結果 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx=0;        
XB=zeros(pl,1); 
for xd=1:xp
    xx=xx+1;
    XB(XP(xd))=MXB(xx);
end

disp('繰返計算回数')
disp(ii)
disp('収束誤差')
disp(Er(ii))
for pp=1:pl
    if XB(pp)<0
    disp('負のパラメータ有り')
    break
    end
end

%%%%%%% t-値の算出 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FH=-inv(HE);
fh=0;
tv=zeros(pl,1);            %--t-値の初期値
for l=1:pl
    if XB(l)~=0
       fh=fh+1;
       tv(l)=XB(l)/double(sqrt(FH(fh,fh)));
    end
end
    
%%%%%%% Beta,t-値の表示 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Beta=zeros(Jmax-1,pk+1);
p=0;
for l=1:Jmax-1
    for M=1:pk+1
        p=p+1;
        Beta(l,M)=XB(p);
    end
end

TV=zeros(Jmax-1,pk+1);
r=0;
for l=1:Jmax-1
    for M=1:pk+1
        r=r+1;
        TV(l,M)=tv(r);
    end
end
disp('ベータの値（段階×変数）')
disp(Beta)
disp('t-値（段階×変数）')
disp(TV)

%%%%%%% ハザード関数,推移確率行列の算出 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ha=zeros(1,Jmax);
for h=1:Jmax-1
    Ha(h)=double(exp(XBP*Beta(h,:)'));
end

Pai=zeros(Jmax,Jmax);
for i=1:Jmax-1
    for j=i:Jmax-1
        if i==j
           Pai(i,i)=Pai(i,i)+double(exp(-Ha(i)*ZP));
        elseif j-i==1
           Pai(i,j)=Pai(i,j)+(Ha(i)/(Ha(i)-Ha(j)))*(-double(exp(-Ha(i)*ZP))+double(exp(-Ha(j)*ZP)));
        else
           P=zeros(1,j-i+1);
           p=0;
           for k=i:j
               Pr=prod(Ha(i:j-1));
               PP=zeros(1,j-i);
               pp=0;
               p=p+1;
               for m=i:j
                   if k~=m
                      pp=pp+1;
                      PP(pp)=PP(pp)+Ha(m)-Ha(k);
                   end
               end
               P(p)=P(p)+(Pr/prod(PP))*double(exp(-Ha(k)*ZP));
           end
           Pai(i,j)=sum(P);
        end
    end
end

for l=1:Jmax
    Pai(l,Jmax)=Pai(l,Jmax)+(1-sum(Pai(l,1:Jmax-1)));
end

disp('ハザード関数')
disp(Ha)
disp('推移確率行列')
disp(Pai)


%%%%%%% 平均レーディング期待寿命 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LI=zeros(1,Jmax-1);
for k=1:K
    XBP=Xk(:,k)';
    Ha=zeros(1,Jmax);
    for h=1:Jmax-1
        Ha(h)=double(exp(XBP*Beta(h,:)'));
    end
    for l=1:Jmax-1
        LI(l)=LI(l)+1/Ha(l);
    end
end
LI=LI/K;

disp('平均レーディング期待寿命')
disp(LI)

%%%%%%% 対数尤度関数,AICの算出 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

L=0;
for k=1:K
    XBP=Xk(:,k)';
    ZP=Ins(k);
        
Ha=zeros(1,Jmax);
for h=1:Jmax-1
    Ha(h)=double(exp(XBP*Beta(h,:)'));
end
if Af(k)==Jmax
Pai=zeros(Jmax,Jmax);
for i=1:Jmax-1
    for j=i:Jmax-1
        if i==j
           Pai(i,i)=Pai(i,i)+double(exp(-Ha(i)*ZP));
        elseif j-i==1
           Pai(i,j)=Pai(i,j)+(Ha(i)/(Ha(i)-Ha(j)))*(-double(exp(-Ha(i)*ZP))+double(exp(-Ha(j)*ZP)));
        else
           P=zeros(1,j-i+1);
           p=0;
           for t=i:j
               Pr=prod(Ha(i:j-1));
               PP=zeros(1,j-i);
               pp=0;
               p=p+1;
               for m=i:j
                   if t~=m
                      pp=pp+1;
                      PP(pp)=PP(pp)+Ha(m)-Ha(t);
                   end
               end
               P(p)=P(p)+(Pr/prod(PP))*double(exp(-Ha(t)*ZP));
           end
           Pai(i,j)=sum(P);
        end
    end
end

for l=1:Jmax
    Pai(l,Jmax)=Pai(l,Jmax)+(1-sum(Pai(l,1:Jmax-1)));
end

L=L+log(Pai(Be(k),Af(k)));
else
   Pai=zeros(Jmax,Jmax);
   i=Be(k);
     j=Af(k);
        if i==j
           Pai(i,i)=Pai(i,i)+double(exp(-Ha(i)*ZP));
        elseif j-i==1
           Pai(i,j)=Pai(i,j)+(Ha(i)/(Ha(i)-Ha(j)))*(-double(exp(-Ha(i)*ZP))+double(exp(-Ha(j)*ZP)));
        else
           P=zeros(1,j-i+1);
           p=0;
           for t=i:j
               Pr=prod(Ha(i:j-1));
               PP=zeros(1,j-i);
               pp=0;
               p=p+1;
               for m=i:j
                   if t~=m
                      pp=pp+1;
                      PP(pp)=PP(pp)+Ha(m)-Ha(t);
                   end
               end
               P(p)=P(p)+(Pr/prod(PP))*double(exp(-Ha(t)*ZP));
           end
           Pai(i,j)=sum(P);
        end
    L=L+log(Pai(Be(k),Af(k)));
end 
end

AIC=-2*L+2*xp;
disp('対数尤度')
disp(L)
disp('AIC')
disp(AIC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
