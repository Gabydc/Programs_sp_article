%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [U,S]=defpodf_Dta(x,dir)
%x=zd;
xm=mean(x,2);

for i=1:size(x,2)
    x(:,i)=x(:,i)-xm;
    if norm(x(:,i))~=0
    x(:,i)=x(:,i)/norm(x(:,i));
    else x(:,i)==0;
    end
        
   
end

lx=size(x,1);
ly=size(x,2);
D=x'*x;
[V,L] = eigs(D,ly);

S=zeros(lx,ly);
g=diag(L);
%g1=sqrt(g.*g);
g1=1./g;
S(1:ly,1:ly)=diag(g1);
mg=max(g1);
S=sparse(S);
V=sparse(V);
x=sparse(x);
U=x*V*S';
%[V,D] = eigs(X,n);
figure(3000)
hd=plot(log(g1/mg),'ob'); 
axis('tight')
title('Eigenvalues X=Z*Z^T','FontSize',16);
ylabel('log(Value) ','FontSize',16)
xlabel('Eigenvalue','FontSize',16)
hold on
figure(4000)
for i=1:lx
hv=plot(V(i,:),'ob'); 
hold on
end
axis('tight')
title('V','-','FontSize',16);
ylabel('log(Value) ','FontSize',16)
xlabel('Eigenvalue','FontSize',16)



file='eig_pod';
B=[dir file '.fig'];
saveas(hd,B)
B=[dir  file '.jpg'];
saveas(hd,B)
file='eig_plot';
B=[dir file '.fig'];
saveas(hv,B)
B=[dir  file '.jpg'];
saveas(hv,B)