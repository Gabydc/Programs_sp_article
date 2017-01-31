
function[xf,ee]=DICCG_te(a,b,xi,iteration,tol,z,l,eig,step,dir,nit,dv,md)
n=length(a);
nn=size(z,2);
e=z'*a*z;
nn=size(e,2);
ei=sparse(inv(e));
if mod(step,dv+1)==0
if eig==1
    if nit==0
        
    if n<3000
     ez=ei*z';
q=z*ez;
pd=eye(n)-a*q;
  pda=pd*a;   
%fprintf('DICCG')
[Va,Da] = eigs(inv(l)*pda*inv(l'),n);
 Da=diag(Da);
 Da=abs(Da);
 [Ve,De] = eigs(e,nn);
 De=diag(De);
 De=abs(De);
  conde=max(De)/min(De);
  figure(90)
 % subplot(2,1,2)
h21=semilogy(De,'or');
axis('tight')

title(['Eigenvalues E, step ' num2str(step) ],'FontSize',16)
ylabel('log(value)','FontSize',16)
xlabel('Eigenvalue','FontSize',16)
file=['eigsE' num2str(step) 'step'];
folder=['eigs'];
mkdir([dir], folder)
dir = [dir folder '/'];
    B=[dir  file  '.fig'];
    saveas(h21,B)
    B=[dir  file  '.jpg'];
    saveas(h21,B)
hold on
 conda=condest(inv(l)*pda*inv(l'),2);
 condadeff=max(Da)/min(Da);
figure(30)
clf
h2=semilogy(Da(1+nn:n),'*');
md
ylim([md(1) max(Da)])
xlim([1+nn n])
title(['Eigenvalues PM^{-1}A, step ' num2str(step)], 'FontSize',16)
ylabel('log(value)','FontSize',16)
xlabel('Eigenvalue','FontSize',16)
hold on
file=['eigsPA' num2str(step) 'step'];
    B=[dir  file  '.fig'];
    saveas(h2,B)
    B=[dir  file  '.jpg'];
    saveas(h2,B)
    end
else h2=0; h21=0;
        end 
end
end
%size(ez)


[pb]=deflatevect(z,ei,a,b);
lb=l\b;
plb=l\pb;
nor1=abs(lb'*lb);
r0=b-a*xi;
[r0]=deflatevect(z,ei,a,r0);
r0=l\r0;
nr0=norm(r0);
p0=l'\r0;
nor=abs(lb'*lb);
%nor=1;
for iter=1:iteration
         [ap]=deflatevect(z,ei,a,a*p0);
         [apt]=tdeflatevect(z,ei,a,p0);
     alpha=(r0'*r0)/((ap)'*p0);   
     xf=xi+alpha*p0;
     r=r0-alpha*(l\ap);
     beta=(r'*r)/(r0'*r0);
     p=l'\r+beta*p0;  
     p0=p;
     r0=r;
     color=[0.3 0.5 0.7];

      ee(iter)=abs(r'*r)/nr0;

       flag=0;

     if (ee(iter)>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     xi=xf;
end
[xf]=tdeflatevect(z,ei,a,xf);
qb=z'*b;
qb=ei*qb;
qb=z*qb;
xf=qb+xf;