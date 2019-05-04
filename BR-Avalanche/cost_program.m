% cost function for multiple objects, multiple failure modes

% number of objects
M=size(obj,1);
% number of mode
N=6;

% read the parameters to matrices
cm=zeros(M,N);
for i=1:N
    v_names{i}=['cm',num2str(i)];
    cm(:,i)=[eval(v_names{i})];
end

alpha=zeros(M,N);
for i=1:N
    v_names{i}=['alpha',num2str(i)];
    alpha(:,i)=[eval(v_names{i})];
end

m=zeros(M,N);
for i=1:N
    v_names{i}=['m',num2str(i)];
    m(:,i)=[eval(v_names{i})];
end

a=zeros(M,N);
for i=1:N
    v_names{i}=['a',num2str(i)];
    a(:,i)=[eval(v_names{i})];
end

clear v_names i


t=0.1:0.1:10;
C=zeros(M,size(t,2));
r=0.02;
ca=[5;5;5];
% run loop through each object
for i=1:M
   
    display(['Object: ',num2str(i)])
    [t,C(i,:)]=cost_function_multimodes(alpha(i,:),m(i,:),r,k0(i),cs3(i),ca(i),cm(i,:),a(i,:),t);
    set(gcf,'name',['object ',num2str(i)]);
    
    
end


        