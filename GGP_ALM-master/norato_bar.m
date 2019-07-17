function [d,dd_dL,dd_dh,dd_dxi,dd_deta]=norato_bar(xi,eta,L,h)
d=((L./2*sqrt(xi.^2./(xi.^2+eta.^2))+sqrt(h.^2/4-L.^2/4*eta.^2./(xi.^2+eta.^2))).*(xi.^2./(xi.^2+eta.^2)>=(L.^2./(h.^2+L.^2)))+h./2.*sqrt(1+xi.^2./eta.^2).*(~(xi.^2./(xi.^2+eta.^2)>=(L.^2./(h.^2+L.^2))))).*(xi~=0|eta~=0)+sqrt(2)/2*h*(xi==0&eta==0);
dd_dL=(1/2*sqrt(xi.^2./(xi.^2+eta.^2))-L*eta.^2/4./(sqrt(h^2/4-L^2/4*eta.^2))./(xi.^2+eta.^2)).*(xi.^2./(xi.^2+eta.^2)>=(L^2/(h^2+L^2))).*(xi~=0&eta~=0);
dd_dh=((h/4./sqrt(h^2/4-L^2/4*eta.^2./(xi.^2+eta.^2))).*(xi.^2./(xi.^2+eta.^2)>=(L^2/(h^2+L^2)))+sqrt(xi.^2./eta.^2+1)/2.*(~(xi.^2./(xi.^2+eta.^2)>=(L^2/(h^2+L^2))))).*(xi~=0&eta~=0);
s1=xi.^2+eta.^2;
dd_dxi=((L/4*(2*xi./(s1)-2*xi.^3./s1.^2)./sqrt(xi.^2./s1)+L^2*eta.^2.*xi./(4*sqrt(h^2/4-L^2/4*eta.^2./s1).*s1.^2)).*(xi.^2./(xi.^2+eta.^2)>=(L^2/(h^2+L^2)))+h*xi./2./eta.^2./sqrt(xi.^2./eta.^2+1).*(~(xi.^2./(xi.^2+eta.^2)>=(L^2/(h^2+L^2))))).*(xi~=0&eta~=0);
dd_deta=(((L^2*eta.^3/2./s1.^2-L^2/2*eta./s1)./(2*sqrt(h^2/4-L^2/4*eta.^2./s1))-L*eta.*xi.^2./(2*sqrt(xi.^2./s1).*s1.^2)).*(xi.^2./(xi.^2+eta.^2)>=(L^2/(h^2+L^2)))-(h*xi.^2./(2*eta.^3.*sqrt(xi.^2./eta.^2+1))).*(~(xi.^2./(xi.^2+eta.^2)>=(L^2/(h^2+L^2))))).*(xi~=0&eta~=0);
