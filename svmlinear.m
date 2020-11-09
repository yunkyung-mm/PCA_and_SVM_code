function [w,bb, acc] = svmlinear_test(dataset, tol,C,maxiter)
%form of dataset is [data variable , response variable]
%data variable is R^n, response variable is only +1 or -1;
%We recommand tol = 10^(-5)
%If it is too small or too large, you can't get result
%
%In linear case, w of output is coefficient x_n
%so we want to get w*x+bb=0
%bb is constant value of solution
%acc is accuracy rate of soultion
%it is 0<=acc<=1 maybe acc value is 1 in separable case
%subject to 0 < alpha < C , default vaule is [] in separable case
%but nonseparable case, input C will control solution
[n, p]= size(dataset);
% sorting dataset according to response variable
% Define matrix a, response variable is +1
% else is matrix b
r=0;
s=0;
for i=1:n;
if dataset(i,3)==1;
 r=r+1;
 a(r,:)=dataset(i,:);
else
 s=s+1;
 b(s,:)=dataset(i,:);
end
end
dataset=[a;b];
%substract response variable
y=dataset(:,p);
%data variable
x=dataset(:,1:p-1);
for i=1:n;
 for j=1:n;
 Q(i,j)=y(i,1)*y(j,1)*(x(i,:)*x(j,:)');
 end
end
34
fun = @(a1) a1'*Q*a1./2-ones(1,n)*a1;
a0=ones(n,1);
y1=[y';-y'; -eye(n)];
toleq=[tol ; tol ; zeros(n,1)];
Aeq=[];
beq=[];
lb=0;
ub=C;
options = optimset('Display','iter','TolFun',10^(-5));
options.MaxIter=maxiter;
a1=fmincon(fun,a0,y1,toleq,Aeq,beq,lb,ub,[],options);
w=(a1.*y)'*x;
bb=((-min(w*x(1:length(a),1:p-1)')+1)+(-max(w*x(length(a)+1:n,1:p1)')-1))/2;
notsep=sum(abs(sign(w*x'+bb)'-y))/2;
acc=1-(notsep/n);
