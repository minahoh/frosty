%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Author: Minah Oh
%Date: 04/27/2020
%Program Description:By using multigrid, this program solves the matrix 
%system Ax=b arising from the Fourier-FEMs applied to the H(div) problem
%on an axisymmetric domain for Fourier modes not equal to zero.
%
%Details can be found in the following paper:
%   https://doi.org/10.1016/j.jmaa.2020.124209
%
%Input: 
% 1. Make mesh1.mat, mesh2.mat, ... , through the final mesh corresponding to
% your domain that contains variables [p,e,t,BdyEdges,new_ele]
%  (1) p,e,t has the same format as those created by using Matlab's "initmesh"
%  (2) BdyEdges is a vector of length the same as the number of edges in the
%      mesh. If the i-th edge is on the boundary then BdyEdges(i)=1 and zero
%      otherwise.
%  (3) new_ele is a matrix of size Number of Elements by 3. The i-th row of
%      new_ele saves the three edge numbers that makes the i-th triangle.
% 2. Save RHS_r.m, RHS_theta.m, RHS_z.m where 
% F(r,z)=(RHS_r(r,z), RHS_theta(r,z), RHS_z(r,z)) is the input function
% 3. n: Fourier_mode of interest not equal to zero.
% 4. mesh_level: The final mesh level to solve the Ax=b system.
%
% Output: The soltuion vector x to Ax=b. (x=main(mesh_level,n))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c2=main(mesh_level, n)



A=cell(mesh_level,1);
b=cell(mesh_level,1);
Mcell=cell(mesh_level,1);
Node=cell(mesh_level,1);
Edge=cell(mesh_level,1);
Ele=cell(mesh_level,1);
New_ele=cell(mesh_level,1);
Basis=cell(mesh_level,1);
Ncell=cell(mesh_level,1);
BdyEdgeCell=cell(mesh_level,1);
VPatchCell1=cell(mesh_level,1);
VPatchCell2=cell(mesh_level,1);

num=7;


for k=1:mesh_level
[Matrix,Vector,Mass,N,E,T,NewT,VP1,VP2,B,Normal,BDYEDGES,n]=Ab(k,n);

A(k)={Matrix};
b(k)={Vector};
Mcell(k)={Mass};
Node(k)={N};
Edge(k)={E};
Ele(k)={T};
New_ele(k)={NewT};
Basis(k)={B};

Ncell(k)={Normal};
BdyEdgeCell(k)={BDYEDGES};

VPatchCell1(k)={VP1};
VPatchCell2(k)={VP2};

end


QQcell=cell(mesh_level,1);

QQcell(1)={0}; %This is for the coarsest grid.


for ii=1:mesh_level-1
    
    normal=cell2mat(Ncell(ii+1));
    
    new_ele=cell2mat(New_ele(ii));
    basis=cell2mat(Basis(ii));
    ele=cell2mat(Ele(ii));
    edge=cell2mat(Edge(ii));
    [N_ele,~]=size(ele);
    [N_edge,~]=size(edge);
    
    node2=cell2mat(Node(ii+1));
    new_ele2=cell2mat(New_ele(ii+1));
    basis2=cell2mat(Basis(ii+1));
    ele2=cell2mat(Ele(ii+1));
    edge2=cell2mat(Edge(ii+1));
    [N_ele2,~]=size(ele2);
    [N_edge2,~]=size(edge2);
    
    BdyEdges2=cell2mat(BdyEdgeCell(ii+1));
  
   prealloc=16*N_ele2;
   i=zeros(prealloc,1);
   j=zeros(prealloc,1);
   s=zeros(prealloc,1);
   
   %new
   prealloc2=10*N_ele2;
   i2=zeros(prealloc2,1);
   j2=zeros(prealloc2,1);
   s2=zeros(prealloc2,1);
    
    
   count=1;
   count2=1;
   for kk2=1:N_ele2
       
        c=mod(kk2,N_ele);
       if c==0
          c=N_ele; 
       end
       
    X1=node2(ele2(kk2,1),1);
    Y1=node2(ele2(kk2,1),2);
    X2=node2(ele2(kk2,2),1);
    Y2=node2(ele2(kk2,2),2);
    X3=node2(ele2(kk2,3),1);
    Y3=node2(ele2(kk2,3),2);
    [X,Y,Wx,Wy]=triquad(num,[X1 Y1;X2 Y2;X3 Y3]);
    

   

    for kk0=1:3 
           for k0=1:3 
    j2(count2)=new_ele(c,kk0);
    i2(count2)=new_ele2(kk2,k0);
   
    
    
    x1=node2(edge2(new_ele2(kk2,k0),1),1);
    y1=node2(edge2(new_ele2(kk2,k0),1),2);
    x2=node2(edge2(new_ele2(kk2,k0),2),1);
    y2=node2(edge2(new_ele2(kk2,k0),2),2);
       
  
 
 
   s2(count2)=(((basis(kk0,1,c)+0.5.*basis(kk0,2,c).*(x1+x2)).*normal(k0,1,kk2)...
               +(basis(kk0,3,c)+0.5.*basis(kk0,2,c).*(y1+y2)).*normal(k0,2,kk2)).*(sqrt((x1-x2).^2+(y1-y2).^2))).*0.5;
           
           
           if  BdyEdges2(new_ele2(kk2,k0))>0 
               s2(count2)=2.*s2(count2);
           end
  
           
           s2(count2)=s2(count2)./(sqrt((x1-x2).^2+(y1-y2).^2));

               
               
           count2=count2+1;
           end
           
    end
           j2(count2)=N_edge+c;
           i2(count2)=N_edge2+kk2;
           s2(count2)=n.*basis(4,2,c); 
           count2=count2+1;
    
    
   
      
       for k0=1:3
           for kk0=1:3
       i(count)=new_ele2(kk2,k0);
       
      
       j(count)=new_ele(c,kk0);
       
       fQ=@(r,z) ((basis2(k0,1,kk2)+basis2(k0,2,kk2).*r).*(basis(kk0,1,c)+basis(kk0,2,c).*r) ...
                + (basis2(k0,1,kk2)+basis2(k0,2,kk2).*r)./n.*(basis(kk0,1,c)+basis(kk0,2,c).*r)./n ...
                + (basis2(k0,3,kk2)+basis2(k0,2,kk2).*z).*(basis(kk0,3,c)+basis(kk0,2,c).*z)).*r;
       
       s(count)=Wx'*feval(fQ,X,Y)*Wy;
       count=count+1;
           end
          i(count)=new_ele2(kk2,k0);
          j(count)=N_edge+c;
          fQ=@(r,z) ((basis2(k0,1,kk2)+basis2(k0,2,kk2).*r)./n.*basis(4,2,c).*r).*r;
          s(count)=Wx'*feval(fQ,X,Y)*Wy;
          count=count+1;
       end
       
       for kk0=1:3
       i(count)=N_edge2+kk2;
       j(count)=new_ele(c,kk0);
       
       fQ=@(r,z) (basis2(4,2,kk2).*r.*(basis(kk0,1,c)+basis(kk0,2,c).*r)./n).*r;
       s(count)=Wx'*feval(fQ,X,Y)*Wy;
       count=count+1;
       end
       
       i(count)=N_edge2+kk2;
       j(count)=N_edge+c;
       fQ=@(r,z) (basis2(4,2,kk2).*r.*basis(4,2,c).*r).*r;
       s(count)=Wx'*feval(fQ,X,Y)*Wy;
       count=count+1;
       
   end
     
   Q=sparse(i2,j2,s2,N_edge2+N_ele2,N_edge+N_ele);
   
   QQcell(ii+1)={Q};
   
   clear i j s q node edge ele new_ele basis node2 edge 2 ele2 new_ele2 basis2
end


 
for ii=mesh_level:mesh_level
    
    V=cell2mat(A(ii));
    vv=cell2mat(b(ii));

    ele=cell2mat(Ele(ii));
    edge=cell2mat(Edge(ii));
    [N_ele,~]=size(ele);
    [N_edge,~]=size(edge);
    
 
        c2=rand(N_edge+N_ele,1);
        u0norm=sqrt(c2'*V*c2);

ee=1;


   while ee>10^-9

   c1=c2;
   c2=mg(ii,c2,vv,n,QQcell,A,VPatchCell1, VPatchCell2, New_ele, Edge);
   er=c1-c2;
   ee=sqrt(er'*V*er)./u0norm;
  
   end

   
end



end
