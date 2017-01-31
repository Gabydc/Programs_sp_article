function [xf,ee,md]=ICCG_te(a,b,xi,iteration,tol,l,dir,eig,step,nit,md)
n=length(a);

if nit==0
if eig==1
    if step==1
%fprintf('ICCG')
%[Va,Da] = eigs(a,n);
[Va,Da] = eigs(inv(l)*a*inv(l'),n);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(a,n);
 condadeff=max(Da)/min(Da);
figure(80)
clf
h2=semilogy(Da,'*');
md(1)=min(Da);
md(2)=max(Da);
axis('tight')
title('Eigenvalues M^{-1}A','FontSize',16)
ylabel('log(value)','FontSize',16)
xlabel('Eigenvalue','FontSize',16)
hold on
h21=0;    
folder=['eigs'];
mkdir([dir], folder)
dir = [dir folder '/'];
file=['eigs' num2str(step) 'step'];
    B=[dir  file  '.fig'];
    saveas(h2,B)
    B=[dir  file  '.jpg'];
    saveas(h2,B)
    %%Eigenvalues of A
    [V,D] = eigs(-a,n);  
figure(81)
clf
h21=semilogy(D,'b*');
m1=min(D);
m2=max(D);
%axis([1 n m1 m2])
axis('tight')
title('Eigenvalues A','FontSize',16)
ylabel('log(value)','FontSize',16)
xlabel('Eigenvalue','FontSize',16)
hold on   
file=['eigs' num2str(step) 'stepA'];
    B=[dir  file  '.fig'];
    saveas(h21(1),B)
    B=[dir  file  '.jpg'];
    saveas(h21(1),B)
    
    
    

else h2=0; h21=0; md=md;
    end 
end
else md=md;
end



r0=b-a*xi;
lb=l\b;
r0=l\r0;
p0=l'\r0;
nor=abs(lb'*lb);
nr0=abs(r0'*r0);
%nor=1;    
text = [dir 'error.txt'];
fileID = fopen(text,'a');
fprintf(fileID,'Iteration & $||M^{-1}b||$ & $||M^{-1}r||$ & $\frac{||M^{-1}r||}{||M^{-1}b||}$ & $max(abs(b-Ax))$ \\\\ \n');
fclose(fileID);
for iter=1:iteration 

    ap=a*p0;
    alpha=(r0'*r0)/((ap)'*p0);
    xf=xi+alpha*p0;
    r=r0-alpha*(l\ap);
    beta=(r'*r)/(r0'*r0);
    p=l'\r+beta*p0;
    p0=p;
    r0=r;
    rr=abs(r'*r);
    ee(iter)=abs(r'*r)/nor;
    ee(iter);
    %pause
    erb=ee(iter);
    e=max(abs(b-a*xf));
    fileID = fopen(text,'a');
    figure(700)
    if nit<5
        subplot(2,2,nit+1)
        hr=plot(step,resNorm,plotStyle{mod(step,6)+1},'Color', colors(mod(step,6)+1,:));
        hold on
        title( ['NR error, NR Iteration ' num2str(nit)] )
    end
    fprintf(fileID,'%d & %1.1e & %1.1e & %1.1e & %1.1e \\\\ \n', iter,nor,rr,erb,e);
    fclose(fileID);
    %pause
    flag=0;
    
    if (ee(iter)>=tol)
        flag=1;
    end
    if flag==0
        break
    end
    xi=xf;

end
end