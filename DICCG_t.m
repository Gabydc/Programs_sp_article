
function[xf,ee]=DICCG_t(a,b,xi,iteration,tol,z,l,eig,step)
n=length(a);
if eig==1
fprintf('ICCG')
[Va,Da] = eigs(a,n);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(a,n)
 condadeff=max(Da)/min(Da)
figure(80)
h2=semilogy(Da,'*');
axis('tight')
title('a')
ylabel('log scale')
xlabel('Eigenvalue')
hold on
h21=0;
file=['eigs' num2str(step) 'step'];
    B=[dir  file  '.fig'];
    saveas(h2,B)
    B=[dir  file  '.jpg'];
    saveas(h2,B)
else h2=0; h21=0;
end

e=z'*a*z;
ei=sparse(inv(e));
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