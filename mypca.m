function [newdata,eig2, k] = mypca (x, perc)
[n,p]=size(x);
x0 = bsxfun(@minus,x,mean(x,1));
%centered
covar=x0'*x0/(n-1); %covariance
[V,D]=eig(covar); %V=eigenvector of covar
D=eig(covar); %D=eigenvalue of covar
eig1=[D';V]; %first row - eigenvalue
%eigenvalue of descending order
[n,p]=size(eig1);
for j=1:p-1;
 for i=1:p-1;
 if eig1(1,i) <= eig1(1,i+1)
 temp=eig1(:,i+1);
 eig1(:,i+1)=eig1(:,i);
 eig1(:,i)=temp;
 end
 end
end
sum(abs(D)); %sum of eigenvalue
%eigenvalue and eigenvector depending on perc
S=0;
k=0;
while S/sum(abs(D)) < perc;
 k=k+1;
 S=S+eig1(1,k);
end
eig2=eig1(:,1:k);
eigweight=eig2(2:end,:);
newdata=x0*eigweight;
