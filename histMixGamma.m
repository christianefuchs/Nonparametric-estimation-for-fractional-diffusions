function []=histMixGamma(T,N,n,a1,b1,a2,b2,u,v,lambda,H,a)
%  --- a=0,1 ---
% a=0 is used to construct the f3 histogram estimators
% a=1 is used to construct the f4 histogram estimators
% u1=gamrnd(a1,b1,[1400 1]); u2=gamrnd(a2,b2,[1400 1]);
% PHII=u*u1+v*u2;
%% Simulating random effects 
rng(231)
for r=1:1
D=T/n; phi=u*gamrnd(a1,b1,[N 1])+v*gamrnd(a2,b2,[N 1]);
W=zeros(n,N); DW=zeros(n,N); X=zeros(n+1,N);
for i=1:N
    [w,t]=fbm1d(H,n,T);
    W(:,i)=w(2:n+1);
end
DW(1,:)=W(1,:);
for  j=2:n
    DW(j,:)=W(j,:)-W(j-1,:);
end
for i=1:N
    for j=2:n+1
        X(j,i)=X(j-1,i)+(-lambda*X(j-1,i)+phi(i))*D+DW(j-1,i);
    end
end
y=(1/T)*X(n+1,:)+(a*lambda*D/T)*sum(X(1:n,:));
%% Histogram construction
x=[0:0.15:7.4]; delta=0.15; p=length(x); f=zeros(p,1); 
for s=1:p-1
    U=0;
    for i=1:N
        if (y(i) >= x(s)) & (y(i) < x(s+1))
        U=U+1;
        end;
    end
    f(s)=U/(N*delta);
end
for s=1:p-1
    plot([x(s); x(s+1); x(s+1)],[f(s); f(s); f(s+1)],'--g')
hold on;
end
end
%% Exact density function
p=length(x2); z=zeros(p,1);
 for s=1:p
 z(s)=u*(1/((b1^a1)*gamma(a1))*x2(s)^(a1-1)*exp(x2(s)/b1))+v*(1/((b2^a2)*gamma(a2))*x2(s)^(a2-1)*exp(x2(s)/b2));
 end
 plot(x,z,'r')



end