%  load variabledf.mat
% character=@(x) Wgp_modified(ugp(:,1),ugp(:,2),x,p,np,nY,Yc);
fun=@(x) compliance(x);
dfun=@(x) dcompliance(x);
% range=[zeros(90,1);ones(90,1)];
range= [1 1 ]';
range=repmat(range,95,1);
size(range);
N=90;
x0= X;
finite_difference_check(fun,dfun,x0,range,N)