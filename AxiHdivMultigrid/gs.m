%%%%%%%%%%%%%%%%%%%%%
%Author: Minah Oh
%Date: 04/27/2020
%%%%%%%%%%%%%%%%%%%%%

%Multiplicative Subspace Correction Method  

function u1=gs(A,b,u0,Vpatch1,Vpatch2) 


[N_node,~]=size(Vpatch1);
N_edge=max(max(Vpatch1));
u1=u0;

for k=1:N_node
    
   au=A*u1; 
    
   a=Vpatch1(k,:);
   a=a(a>0);
   
   c=Vpatch2(k,:);
   c=c(c>0);
   c=N_edge+c;
   
   d=[a,c];
   
   B=b(d); 
   AU=au(d); 
   SubA=A(d,d); 
   
   diff=B-AU;
   invAdiff=SubA\diff;
   u1(d)=u1(d)+invAdiff;
   
   clear B AU SubA diff invAdiff a c d au
end
%{
for k=1:N

s=A(k,:)*u1;    
u1(k)=u1(k)+(b(k)-s)./A(k,k);

end
%}
