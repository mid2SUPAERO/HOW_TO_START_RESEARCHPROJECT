%  load variabledf.mat
% character=@(x) Wgp_modified(ugp(:,1),ugp(:,2),x,p,np,nY,Yc);
% fun=@(x) compliance(x);
% dfun=@(x) dcompliance(x);
fun=@(x) myfunc(x,p);
dfun=@(x) myder(x,p);
% range=[zeros(90,1);ones(90,1)];
range= [1]';
range=repmat(range,length(TV),1);
size(range);
N=20;
x0= TV;
finite_difference_check(fun,dfun,x0,range,N)