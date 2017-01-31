A =[3 2 2;
    2 3 -2]
A=A';
AAT=A*A';
ATA=A'*A;
[V,L]=eigs(ATA)
S=zeros(3,3);
g=diag(L)
g1=1./g1
S(1:2,1:2)=diag(g1)