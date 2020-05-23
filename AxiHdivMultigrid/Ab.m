%%%%%%%%%%%%%%%%%%%%%
%Author: Minah Oh
%Date: 04/27/2020
%%%%%%%%%%%%%%%%%%%%%

function [A,b,M,node,edge,ele,new_ele,Vpatch1,Vpatch2,basis,normal,BdyEdges,n]=Ab(mesh_level,n)

load(['mesh',num2str(mesh_level),'.mat']);    

num_quadrature=7; 

[~,N_node]=size(p); 
node=p'; 
[~,N_ele]=size(t); 
ele=t(1:3,1:N_ele);
ele=ele'; 
TR=triangulation(ele,node); 
edge=edges(TR); 
[N_edge,~]=size(edge);  
  
%Vpatch1 contains the edge numbers connected to each vertex.
Vpatch1=-ones(N_node,6);
Count4Vpatch=ones(N_node,1);
for k=1:N_edge
    Vpatch1(edge(k,1),Count4Vpatch(edge(k,1)))=k;
    Count4Vpatch(edge(k,1))=Count4Vpatch(edge(k,1))+1;
    
    Vpatch1(edge(k,2),Count4Vpatch(edge(k,2)))=k;
    Count4Vpatch(edge(k,2))=Count4Vpatch(edge(k,2))+1;
   
end

%Vpatch2 contains the element numbers connected to each vertex.
Vpatch2=-ones(N_node,6);
Count4Vpatch=ones(N_node,1);
for k=1:N_ele
    Vpatch2(ele(k,1),Count4Vpatch(ele(k,1)))=k;
    Count4Vpatch(ele(k,1))=Count4Vpatch(ele(k,1))+1;
    
    Vpatch2(ele(k,2),Count4Vpatch(ele(k,2)))=k;
    Count4Vpatch(ele(k,2))=Count4Vpatch(ele(k,2))+1;
    
    Vpatch2(ele(k,3),Count4Vpatch(ele(k,3)))=k;
    Count4Vpatch(ele(k,3))=Count4Vpatch(ele(k,3))+1;
   
end

prealloc=N_ele*16;
prealloc2=N_ele*4;
i=zeros(1,prealloc);
j=zeros(1,prealloc);
s=zeros(1,prealloc);
s_mass=zeros(1,prealloc);
i_rhs=zeros(1,prealloc2);
j_rhs=zeros(1,prealloc2);
s_rhs=zeros(1,prealloc2);
k_rhs=1;

[basis,normal]=rDIVbasis(p,t,new_ele);

k_global=1; 

for k=1:N_ele
    
    X1=node(ele(k,1),1);
    Y1=node(ele(k,1),2);
    X2=node(ele(k,2),1);
    Y2=node(ele(k,2),2);
    X3=node(ele(k,3),1);
    Y3=node(ele(k,3),2);
    
   [X,Y,Wx,Wy]=triquad(num_quadrature,[X1 Y1;X2 Y2;X3 Y3]);

   for k0=1:3
      
      f=@(r,z) ((basis(k0,1,k)+basis(k0,2,k).*r).^2+(basis(k0,3,k)+basis(k0,2,k).*z).^2+((basis(k0,1,k)+basis(k0,2,k).*r)./n).^2+(2.*basis(k0,2,k)).^2).*r; 
      s(k_global)=Wx'*feval(f,X,Y)*Wy;  
      
      f_div=@(r,z) ((2.*basis(k0,2,k)).^2).*r;
      s_mass(k_global)=s(k_global)-Wx'*feval(f_div,X,Y)*Wy;
      
      i(k_global)=new_ele(k,k0); 
      j(k_global)=new_ele(k,k0); 
      k_global=k_global+1;
      
      for k1=k0+1:3
      f=@(r,z) ((basis(k0,1,k)+basis(k0,2,k).*r).*(basis(k1,1,k)+basis(k1,2,k).*r)+(basis(k0,3,k)+basis(k0,2,k).*z).*(basis(k1,3,k)+basis(k1,2,k).*z)+((basis(k0,1,k)+basis(k0,2,k).*r)./n).*((basis(k1,1,k)+basis(k1,2,k).*r)./n)+(4.*basis(k0,2,k).*basis(k1,2,k))).*r; 
      s(k_global)=Wx'*feval(f,X,Y)*Wy;  
      
      f_div=@(r,z) (4.*basis(k0,2,k).*basis(k1,2,k)).*r;
      s_mass(k_global)=s(k_global)-Wx'*feval(f_div,X,Y)*Wy; 
      
      i(k_global)=new_ele(k,k0); 
      j(k_global)=new_ele(k,k1); 
      k_global=k_global+1;

      i(k_global)=new_ele(k,k1); 
      j(k_global)=new_ele(k,k0);
      s(k_global)=s(k_global-1);
      s_mass(k_global)=s_mass(k_global-1);
      k_global=k_global+1;
      end
      
      i(k_global)=k+N_edge;
      j(k_global)=new_ele(k,k0);
      f=@(r,z) (((basis(k0,1,k)+basis(k0,2,k).*r)./n).*(basis(4,2,k).*r)+(2.*basis(k0,2,k).*(-n.*basis(4,2,k)))).*r;
      s(k_global)=Wx'*feval(f,X,Y)*Wy;
      
      f_div=@(r,z) (2.*basis(k0,2,k).*(-n.*basis(4,2,k))).*r;
      s_mass(k_global)=s(k_global)-Wx'*feval(f_div,X,Y)*Wy; 
      
      k_global=k_global+1;

      j(k_global)=k+N_edge;
      i(k_global)=new_ele(k,k0);
      s(k_global)=s(k_global-1);
      s_mass(k_global)=s_mass(k_global-1);
      k_global=k_global+1;
   end
   
   i(k_global)=k+N_edge;
   j(k_global)=k+N_edge;
   f=@(r,z) ((basis(4,2,k).*r).^2+n.^2.*basis(4,2,k).^2).*r;
   s(k_global)=Wx'*feval(f,X,Y)*Wy;
   
   f_div=@(r,z) (n.^2.*basis(4,2,k).^2).*r;
   s_mass(k_global)=s(k_global)-Wx'*feval(f_div,X,Y)*Wy; 
   
   k_global=k_global+1;

   for k0=1:3
   i_rhs(k_rhs)=new_ele(k,k0);
   j_rhs(k_rhs)=1;
   f_rhs=@(r,z) ((basis(k0,1,k)+basis(k0,2,k).*r).*RHS_r(r,z)+(basis(k0,3,k)+basis(k0,2,k).*z).*RHS_z(r,z)+(basis(k0,1,k)+basis(k0,2,k).*r)./n.*RHS_theta(r,z)).*r;
   s_rhs(k_rhs)=Wx'*feval(f_rhs,X,Y)*Wy;
   k_rhs=k_rhs+1;
   end
   i_rhs(k_rhs)=k+N_edge;
   j_rhs(k_rhs)=1;
   f_rhs=@(r,z) (basis(4,2,k).*r.*RHS_theta(r,z)).*r;
   s_rhs(k_rhs)=Wx'*feval(f_rhs,X,Y)*Wy;
   k_rhs=k_rhs+1;
end


A=sparse(i,j,s,N_edge+N_ele, N_edge+N_ele);
b=sparse(i_rhs,j_rhs,s_rhs,N_edge+N_ele,1);
M=sparse(i,j,s_mass,N_edge+N_ele, N_edge+N_ele);

end

