% solve the equation q(T)=0
% example:
% 
% t=0.1:0.1:30;
% [t,f,C,result]=solveqT(0.05,1.5,0.02,0.2,10,100,50,t);
% 
% the results will be:
% 
% the optimal T=10.6
% the cost C=913.0527
function [t,f,C,result]=solveqT(a,m,r,K,Cs,Ca,Cm,t)

% a is alpha
% m is m
% r is rho
% K is Ko
% Cs, Ca, Cm are the costs
% t is the considered period of time, e.g. t=0:0.1:30
% C is the cost function
% f is the derivative of the cost function (f=q(T))

% memory pre-allocation
f=zeros(1,size(t,2));
C=zeros(1,size(t,2));
% calculation of f and C as function of time T
for i=1:size(t,2)
    
    f1=@(x) (((m*a.*(x.^(m-1)).*exp(-a.*(x.^m)))).^2);
    v1=integral(f1,0,inf);
    f2=@(x) (exp(a.*(x.^m))-r.*x);
    v2=integral(f2,0,t(i));
    f(i)= v1*Cm*exp(a*(t(i)^m))*(1-exp(-r*t(i))) - Cm*r*v2 *v1 - r* (Cs+Ca);
    C(i)=(exp(-r*t(i))*(Cs+Ca)+ K*(1-exp(-r*t(i)))/r +Cm*v2*v1 )/(1-exp(-r*t(i)));
end

% find the last instant t that the derivative f(t) is negative, i.e. the
% cost funtion tends to decrease.
ind=find(f<0);
j=ind(end);

display(['the optimal T=',num2str(t(j))])

display(['the cost C=',num2str(C(j))])

figure();

% PropertiesFigures;
plot(t,f);
xlabel('T')
ylabel('q(T)');

grid on; box on;

figure();

% PropertiesFigures;
plot(t,C);
xlabel('T')
ylabel('cost function C(T)');
legend('a Nam')
grid on; box on;

% final result: 3 columns [time, derivative,cost]
result=[t',f',C'];

end
