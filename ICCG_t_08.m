function [xf,ex]=ICCG_t_08(a,b,xi,iteration,tol,l,dir)
%a=A;

n=length(a);
xii=xi;
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
ex(iter)=norm(xf-xi)/norm(xi);
     figure(500)
     color=[0.2 0.8 0.6];
     h=semilogy(iter,ex(iter),'o','Color',color);
     hold on
ylabel('log(||M^{-1}r^k||_2/||M^{-1}b||_2)','FontSize',16)
xlabel('Iteration','FontSize',16)
fprintf(fileID,'%d & %1.1e & %1.1e & %1.1e & %1.1e \\\\ \n', iter,nor,rr,erb,e);
fclose(fileID);  
    %pause
    flag=0;
    
    if (ex(iter)>=tol)
        flag=1;
    end
    if flag==0
        break
    end
    xi=xf;

end

end
