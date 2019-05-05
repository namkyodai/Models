% clear all
% 
% %%%%%%% ���v�ʂ̒��o %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Beta=dlmread('alphasamp2.dat');
% % DB=dlmread('hiddenpo.txt');
% ZP=2;
% DB=xlsread('�k���ڎ��f�[�^.xlsx','�k���ڎ�DB','a2:i2536')';    %--�t�@�C���̓ǂݍ���
% Cd=[];                     %--�����ϐ��Ƃ��č̗p����ϐ�
% if length(Cd)==0
% pk=0;                                               %�̂̒�
% [TPP,K]=size(DB);                                   %�łݐ�
% Xk=ones(1,K);                                       %�@�ȍ�
% else
% pk=length(Cd);               %--�萔�ȊO�̐����ϐ��̌�
% for i=1:pk
%      X(i,:)=DB(Cd(i),:)/max(DB(Cd(i),:));     % Explanation Value
% end
% 
% [M,K]=size(X);               %-- �T���v���T�C�Y�FK�̐����グ
% X1=ones(1,K);                %-- �萔���̐ݒ�
% Xk=[X1;X];                   %-- �萔���{���̑��̕ϐ�  
% end
% % %%--1=�ÎR�C2=�O���C3=���c�C4=�R���C5=���R
% % %%--6=���R�C7=�L���C8=���R�C9=�Ďq�C10=�l�c�C11=���]
% % In=6;
% % MNo=0;
% 
% Jmax=3;                       %--�i�K��
% 
% Be=DB(6,:);                  % Before
% Af=DB(8,:);                  % After
% Ins=DB(9,:);                 % Interval
% 
% % Zi=30;                          %-- �̗p�N��
% Ni=15;                            %-- ���e��
% Ti=30;                           %-- �O��_���܂ł̌o�ߔN��
% Xp=[1 1 1];                 %-- �̗p�ϐ�
Beta=samP;
[PPOINT,i]=size(Beta);
% Cd=[7 8];
xp=i;
IBURN=1000;                          %-- �o�[���C����
NPOINT=PPOINT-IBURN;                        %-- �W�{��
% PPOINT=NPOINT+IBURN;
% Bt=DB(:,2)';                  % ���p�J�n�N
% At=DB(:,3)';                  % �����_
% Et=At-Bt;                  % �_���Ԋu
% Et0=DB(:,5)';                  % ���p����
% Id=DB(:,7)';                  % �����̑�����
% 
% pk=length(Cd);
% if pk==0
%     PM1=0;
%     PK=length(Ets);
%     PX1=ones(1,PK);
%     rPXk=[PX1];
% else
% for i=1:pk
%      PX(i,:)=DB(:,Cd(i))';     % Explanation Value
% end
% [PM1,PK]=size(PX);         %-- �T���v���T�C�Y�FK�̐����グ
% PX1=ones(1,PK);          %-- �萔���̐ݒ�
% rPXk=[PX1;PX];             %-- �萔���{���̑��̃x�[�^
% end
% for i=1:PM1+1
%      PXk(i,:)=rPXk(i,:)/max(rPXk(i,:));
% end
% 


Beta0=Beta(IBURN+1:PPOINT,:);
plot(Beta(:,1:i))
% figure(2),plot(Beta0(:,1:i))
[a b]=size(Beta);

for pal=1:b
    Beta1(:,pal)=sort(Beta(:,pal));
end
Beta0=sort(Beta0);
B050=mean(Beta0);
B005=Beta0(round(0.05*NPOINT),:);
B095=Beta0(round(0.95*NPOINT),:);
BB=[B050;B005;B095];


%%%%%%% geweke����� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ge=zeros(1,b);

for m=1:b
    vvar=PPOINT;
    vund=IBURN;
    v1=0.1*(vvar-vund);
    v2=0.5*(vvar-vund);
    
    be1=mean(Beta(vund+1:vund+v1,m));
    be2=mean(Beta(vvar-v2+1:vvar,m));
    ome1=var(Beta(vund+1:vund+v1,m));
    ome2=var(Beta(vvar-v2+1:vvar,m));
    SS1=0;
    for s=1:20
        SV1=0;
        for v=vund+s+1:vund+v1
            SV1=SV1+(Beta(v,m)-be1)*(Beta(v-s,m)-be1);
        end
        SS1=SS1+(1-(s/21))*SV1;
    end
    F1=(ome1+2*SS1)/v1;
    
    SS2=0;
    for s=1:20
        SV2=0;
        for v=vvar-v2+s+1:vvar
            SV2=SV2+(Beta(v,m)-be2)*(Beta(v-s,m)-be2);
        end
        SS2=SS2+(1-(s/21))*SV2;
    end
    F2=(ome2+2*SS2)/v2;
    
    ge(m)=(be1-be2)/sqrt(F1+F2);
end






%%%%%%%%%% �ΐ��ޓx�̎Z�o %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('�x�[�^ ���ρ�')
disp([BB(1,:)])
disp('�x�[�^ 95%�M����ԁ�')
disp([BB(2:3,:)])
disp('geweke����ʁ�')
disp(ge)

%%%%%%% �ΐ��ޓx�֐�,AIC�̎Z�o %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% 
% L=0;
% for k=1:K
%     XBP=Xk(:,k)';
%     ZP=Ins(k);
%         
% Ha=zeros(1,Jmax);
% for h=1:Jmax-1
%     Ha(h)=double(exp(XBP*B050((pk+1)*(h-1)+1:(pk+1)*h)'));
% end
% if Af(k)==Jmax
% Pai=zeros(Jmax,Jmax);
% for i=1:Jmax-1
%     for j=i:Jmax-1
%         if i==j
%            Pai(i,i)=Pai(i,i)+double(exp(-Ha(i)*ZP));
%         elseif j-i==1
%            Pai(i,j)=Pai(i,j)+(Ha(i)/(Ha(i)-Ha(j)))*(-double(exp(-Ha(i)*ZP))+double(exp(-Ha(j)*ZP)));
%         else
%            P=zeros(1,j-i+1);
%            p=0;
%            for t=i:j
%                Pr=prod(Ha(i:j-1));
%                PP=zeros(1,j-i);
%                pp=0;
%                p=p+1;
%                for m=i:j
%                    if t~=m
%                       pp=pp+1;
%                       PP(pp)=PP(pp)+Ha(m)-Ha(t);
%                    end
%                end
%                P(p)=P(p)+(Pr/prod(PP))*double(exp(-Ha(t)*ZP));
%            end
%            Pai(i,j)=sum(P);
%         end
%     end
% end
% 
% for l=1:Jmax
%     Pai(l,Jmax)=Pai(l,Jmax)+(1-sum(Pai(l,1:Jmax-1)));
% end
% 
% L=L+log(Pai(Be(k),Af(k)));
% else
%    Pai=zeros(Jmax,Jmax);
%    i=Be(k);
%      j=Af(k);
%         if i==j
%            Pai(i,i)=Pai(i,i)+double(exp(-Ha(i)*ZP));
%         elseif j-i==1
%            Pai(i,j)=Pai(i,j)+(Ha(i)/(Ha(i)-Ha(j)))*(-double(exp(-Ha(i)*ZP))+double(exp(-Ha(j)*ZP)));
%         else
%            P=zeros(1,j-i+1);
%            p=0;
%            for t=i:j
%                Pr=prod(Ha(i:j-1));
%                PP=zeros(1,j-i);
%                pp=0;
%                p=p+1;
%                for m=i:j
%                    if t~=m
%                       pp=pp+1;
%                       PP(pp)=PP(pp)+Ha(m)-Ha(t);
%                    end
%                end
%                P(p)=P(p)+(Pr/prod(PP))*double(exp(-Ha(t)*ZP));
%            end
%            Pai(i,j)=sum(P);
%         end
%     L=L+log(Pai(Be(k),Af(k)));
% end 
% end
% 
% AIC=-2*L+2*xp;
% disp('�ΐ��ޓx')
% disp(L)
% disp('AIC')
% disp(AIC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% j=zeros(a,Jmax-1);
% for c=1:a
%     for d=1:Jmax-1
%     j(c,d)=1/exp(Beta(c,(d-1)*(pk+1)+1)+Beta(c,(d-1)*(pk+1)+2)*0);
%     if d==2
%         j(c,d)=j(c,1)+j(c,d);
%     end
%     end
% end
% 
% 
% j0=sort(j);
% j050=mean(j0);
% j005=j0(round(0.05*NPOINT),:);
% j095=j0(round(0.95*NPOINT),:);
% jj=[j050;j005;j095];




