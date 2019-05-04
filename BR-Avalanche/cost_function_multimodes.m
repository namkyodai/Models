
% % calculate the cost of ONE object with MULTIPLE failure modes

function [t,C]=cost_function_multimodes(alpha,m,r,K,Cs,Ca,Cm,a,t)
% alpha
% a
% m
% r is rho
% K is Ko
% Cs, Ca, Cm are the costs

% % alpha, m, a, Cm, alpha are ROW vectors

% t is the considered period of time, e.g. t=0:0.1:10
% C is the cost function
% f is the derivative of the cost function (f=q(T))

% memory pre-allocation
f=zeros(1,size(t,2));
C=zeros(1,size(t,2));
% calculation of f and C as function of time T
for i=1:size(t,2)
    for j=1:size(alpha,2)
        f1=@(x) (((m(j)*alpha(j).*(x.^(m(j)-1)).*exp(-alpha(j).*(x.^m(j))))).^2);
        v1(j)=integral(f1,0,inf);
        f2=@(x) (exp(alpha(j).*(x.^m(j)))-r.*x);
        v2(j)=integral(f2,0,t(i));
    end

%     f(i)= v1*Cm*exp(a*(t(i)^m))*(1-exp(-r*t(i))) - Cm*r*v2 *v1 - r* (Cs+Ca);
    C(i)=(exp(-r*t(i))*(Cs+Ca)+ K*(1-exp(-r*t(i)))/r +sum(a.*Cm.*v2.*v1) )/(1-exp(-r*t(i)));
end

% find minimal cost C(j)
j=find(C==min(C));

display(['the optimal T=',num2str(t(j))])

display(['the minimal cost C=',num2str(C(j))])

% plot the cost function
figure();

%PropertiesFigures;
h=plot(t,C);
set(h,'linewidth',2,'color','b');
xlabel('T')
ylabel('cost function C(T)');
% xlim([0,30]);
% ylim([0,40000]);
legend('cost')

grid on; box on;

end
