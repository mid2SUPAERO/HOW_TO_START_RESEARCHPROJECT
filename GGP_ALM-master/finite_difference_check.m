function [error]= finite_difference_check(f,df,x0,range,N)
% Finite Difference Check
% INPUT: 
% f:function handle of the function
% df: function handle of the sensitivity
% x0: point where the sensitivity is computed
% range: range of input variable to scale perturbations
% N: number of component of sensitivity to be tested
N=min(N,length(x0));
% Range_test=randperm(length(x0),N);
df0=df(x0);
Range_test=find(df0);
exponent=-[3,4,5,6];
error=zeros(length(exponent),length(Range_test));
l=0;
for ex=exponent
    l=l+1;
    parfor sk=1:length(Range_test)
        k=Range_test(sk);
        pert=zeros(size(x0));
        pert(k)=range(k)*10.^ex;
        xp=x0+pert/2;
        fp=f(xp);
        xm=x0-pert/2;
        fm=f(xm);
        df_finite_difference=(fp-fm)/(range(k)*10.^ex);
        error(l,sk)=norm(df_finite_difference-df0(k),1)/max(norm(df0(k)),1e-2)*100;
    end
    legend_name{l}=['perturbation = range \times 10^{' ,num2str(exponent(l)),'}'];
end
bar(Range_test,error');
b=gca;
b.YScale = 'log';
grid on
legend(legend_name,'Location','bestoutside')
title('Relative error [%]')
xlabel('sensitivity component')
