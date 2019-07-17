function [dalmc] = dalm_test(X)
Xc=X;
Xk=Xc(1:2:end-2*np);        % Extraire tous les Xk
Lk=Xc(2:2:end-2*np);        % Extraire tous les Lk
Xk=reshape(Xk,nY,np);
Lk=reshape(Lk,nY,np);
% xk & dxk_dXk
xk=Xk(1:end-1,:);
dxk_dXk=zeros(size(xk,1),size(Xk,1));
dxk_dXk(1:end,1:(end-1))=eye(size(xk,1));
dxk_dXk=sparse(dxk_dXk);
% xk+1 & dxk+1_dXk
xk1=Xk(2:end,:);
dxk1_dXk=zeros(size(xk1,1),size(Xk,1));
dxk1_dXk(1:end,2:(end))=eye(size(xk1,1));
dxk1_dXk=sparse(dxk1_dXk);
% lk & dlk_dLk
lk=Lk(1:end-1,:);
dlk_dLk=zeros(size(lk,1),size(Lk,1));
dlk_dLk(1:end,1:(end-1))=eye(size(lk,1));
dlk_dLk=sparse(dlk_dLk);
% lk+1 & dlk+1_dLk
lk1=Lk(2:end,:);
dlk1_dLk=zeros(size(lk,1),size(Lk,1));
dlk1_dLk(1:end,2:(end))=eye(size(lk1,1));
dlk1_dLk=sparse(dlk1_dLk);


uuk=xk+lk/2;
uuk1=xk1+lk1/2;
llk=xk-lk/2;
llk1=xk1-lk1/2;
tanalpha= ((uuk1 - uuk)/nely)*(nY-1);
tanbeta= ((llk - llk1)/nely)*(nY-1);
tv1 = (tanalpha(:) - tand(35))/tand(35)*100;
tv2= (tanbeta(:) - tand(35))/tand(35)*100;
TV=[tv1;tv2];
[almc,dalmc_dtv]=Aggregation_Pi(TV,p);
dtv1_dXk = (100/(tand(35)))*(((nY-1)/nely)* (dxk1_dXk - dxk_dXk));
dtv1_dXk=repmat(dtv1_dXk,np,1);
j1=reshape(repmat(1:np,nY-1,1),[],1);
j2=nY*(repmat(j1,1,nY)-1)+repmat(1:nY,length(j1),1);
i2=repmat((1:(np*(nY-1)))',1,nY);
dtv1_dXk=sparse(i2(:),j2(:),dtv1_dXk(:),(nY-1)*np,nY*np);

% dalmc_dtv1_dXk= reshape(dalmc_dtv1,[],np)'* dtv1_dXk ;

dtv1_dLk = (100/(tand(35)))*(((nY-1)/nely)* (-dlk1_dLk/2 - dlk_dLk/2));
dtv1_dLk=repmat(dtv1_dLk,np,1);
dtv1_dLk=sparse(i2(:),j2(:),dtv1_dLk(:),(nY-1)*np,nY*np);




dtv2_dXk = (100/(tand(35)))*(((nY-1)/nely)* (dxk_dXk - dxk1_dXk));
dtv2_dXk=repmat(dtv2_dXk,np,1);
dtv2_dXk=sparse(i2(:),j2(:),dtv2_dXk(:),(nY-1)*np,nY*np);


dtv2_dLk = (100/(tand(35)))*(((nY-1)/nely)* (-dlk_dLk/2 + dlk1_dLk/2));
dtv2_dLk=repmat(dtv2_dLk,np,1);
dtv2_dLk=sparse(i2(:),j2(:),dtv2_dLk(:),(nY-1)*np,nY*np);
dtv_dXk=[dtv1_dXk;dtv2_dXk];
dtv_dLk=[dtv1_dLk;dtv2_dLk];
dalmc_dXk=dalmc_dtv(:)'*dtv_dXk;
dalmc_dLk=dalmc_dtv(:)'*dtv_dLk;
dalmc=zeros(size(Xc));
dalmc(1:2:end-2*np)=dalmc_dXk;
dalmc(2:2:end-2*np)=dalmc_dLk;

end

