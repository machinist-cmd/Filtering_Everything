function [xk_a,Pk_a,A] = fun_UKF(xk_a,Pk_a,Fk,Gk,Z_true,Qk,sigma_r,sigma_b,sigma_e, xp, nx,nz)
%UKF Filter
%%
zk=Z_true;% 雷达量测:rm bm em
n=nx;%状态维数，CV6,CT为6 CA为9
nz=nz;% 测量维数
Pk=Pk_a([1:nx],[1:nx]);
xk=xk_a([1:nx]);
Pk;
% 提出奇异Pk
d=eig(Pk);
if all(d > 0) == 1
    Pk=Pk;
else 
    Pk=diag(diag(Pk)) ;
end
% xkk=Fk*xk;
% Pkk=Fk*Pk*Fk'+Gk*Qk*Gk';
%UT transformation
alpha=0.98;
kk=0;
beta=2; 
lambda=alpha^2*(n+kk)-n;
Wm=[lambda/(lambda+n),  (0.5/(lambda+n))+zeros(1,2*n)];%权值确定
Wc=[lambda/(lambda+n)+1-alpha^2+beta,   (0.5/(lambda+n))+zeros(1,2*n)];
%产生xk的Sigma点
SPk=sqrt(n+lambda)*(chol(Pk))';
Xsigma0=xk;
for i=1:n 
    Xsigma1(:,i)=xk+SPk(:,i);
    Xsigma2(:,i)=xk-SPk(:,i); 
end
Xsigma=[Xsigma0,Xsigma1,Xsigma2];
%产生xkk的Gama点
for i=1:2*n+1
    Xgama(:,i)=Fk*Xsigma(:,i);
end
xkk=Xgama*Wm';
Pkk=zeros(nx,nx);for i=1:2*n+1; Pkk=Pkk+Wc(i)*((Xgama(:,i)-xkk)*(Xgama(:,i)-xkk)');end; 
Pkk=Pkk+Gk*Qk*Gk';

%产生xkk的Sigma点
SPkk=sqrt(n+lambda)*(chol(Pkk))';
Zsigma0=xkk;
for i=1:n 
    Zsigma1(:,i)=xkk+SPkk(:,i);
    Zsigma2(:,i)=xkk-SPkk(:,i);
end
Zsigma=[Zsigma0,Zsigma1,Zsigma2];
% zkk
for i=1:2*n+1
    [z1,z2,z3] = measurements(Zsigma(:,i), xp); Zgama(:,i)=[z1,z2,z3]';
end

zkk=Zgama*Wm';

Rk=diag([sigma_r^2, sigma_b^2, sigma_e^2]);
Sk=zeros(nz,nz);for i=1:2*n+1; Sk=Sk+Wc(i)*((Zgama(:,i)-zkk)*(Zgama(:,i)-zkk)');end; 
Sk=Sk+Rk;
Ck=zeros(nx,nz);for i=1:2*n+1; Ck=Ck+Wc(i)*((Xgama(:,i)-xkk)*(Zgama(:,i)-zkk)');end;
% 更新
Kk=Ck*inv(Sk);
xk=xkk+Kk* (zk-zkk) ;
Pk=Pkk-Kk*(Sk)*Kk';
A=(2*pi)^(-1/2)*(det(Sk))^(-1/2)*exp((-1/2)*(zk-zkk)'*inv(Sk)*(zk-zkk));


Pk_a([1:nx],[1:nx])=Pk;
xk_a([1:nx])=xk;

end