%%%%%%%%%%%%%%%%%%%%%
%Author: Minah Oh
%Date: 04/27/2020
%%%%%%%%%%%%%%%%%%%%%

function [basis,normal]=rDIVbasis(p,t,new_ele)

node=p';  
[~,N_ele]=size(t); 
ele=t(1:3,1:N_ele);
ele=ele'; 
TR=triangulation(ele,node); 
edge=edges(TR); 
[N_edge,~]=size(edge); 

basis=zeros(4,3,N_ele);
normal=zeros(3,2,N_ele);
abc=zeros(3,3,N_ele);

edge_length=zeros(N_edge,1);
for k_edge=1:N_edge
   edge_length(k_edge)=sqrt((node(edge(k_edge,1),1)-node(edge(k_edge,2),1))^2+(node(edge(k_edge,1),2)-node(edge(k_edge,2),2))^2); 
end


for k=1:N_ele
    
    X1=node(ele(k,1),1);
    Y1=node(ele(k,1),2);
    X2=node(ele(k,2),1);
    Y2=node(ele(k,2),2);
    X3=node(ele(k,3),1);
    Y3=node(ele(k,3),2);
    
    area=polyarea([X1;X2;X3],[Y1;Y2;Y3]);
    
    S=[X1,Y1,1; X2, Y2, 1; X3,Y3,1]; 
       
   %abc is the matrix that saves the local lambda coefficients.
   abc(1,:,k)=S\[1;0;0];
   abc(2,:,k)=S\[0;1;0];
   abc(3,:,k)=S\[0;0;1];
   

    %First Edge Basis 
basis(1,2,k)=edge_length(new_ele(k,1))./(2.*area);
basis(1,1,k)=-basis(1,2,k).*X3;
basis(1,3,k)=-basis(1,2,k).*Y3; 
 
    
     
    midpoint=[(X1+X2)*0.5;(Y1+Y2)*0.5];
   
   
    Normal=zeros(2,1);
    
    Normal(1)=Y1-Y2;
    Normal(2)=X2-X1;
    Nlength=norm(Normal);
    Normal=Normal./Nlength;
    

    
    
    if distance2(midpoint+Normal,[X3;Y3])<distance2(midpoint-Normal,[X3;Y3])
        Normal=-Normal;
     
    end
 

  
  if Normal(1)==0
     if Normal(2)<0
        basis(1,:,k)=-basis(1,:,k); 
       
     end
  end
    
  if Normal(2)==0
     if Normal(1)<0
        basis(1,:,k)=-basis(1,:,k); 
       
     end
  end
  
  if abs(Normal(1))==abs(Normal(2))
     if Normal(1)<0
        basis(1,:,k)=-basis(1,:,k); 
        
     end
  end 
  Normal2=Normal;
  
  %%Normal2 for prolongation
    if Normal(1)==0
     if Normal(2)<0
       Normal2=-Normal;
       
     end
    end
    
  if Normal(2)==0
     if Normal(1)<0
       Normal2=-Normal;
       
     end
  end
  
  if abs(Normal(1))==abs(Normal(2))
     if Normal(1)<0
         
        Normal2=-Normal;
        
     end
  end 

normal(1,1,k)=Normal2(1);
normal(1,2,k)=Normal2(2);
  
  %Second Edge Basis

 basis(2,2,k)=edge_length(new_ele(k,2))/(2*area);
 basis(2,1,k)=-basis(2,2,k)*X1;
 basis(2,3,k)=-basis(2,2,k)*Y1;
      
  
  
    midpoint=[(X2+X3)*0.5;(Y2+Y3)*0.5];
  
        
    Normal(1)=Y2-Y3;
    Normal(2)=X3-X2;
    Nlength=norm(Normal);
    Normal=Normal./Nlength;
    
    
    if distance2(midpoint+Normal,[X2;Y2])<distance2(midpoint-Normal,[X1;Y1])
        Normal=-Normal;
    end
    
  if Normal(1)==0
     if Normal(2)<0
        basis(2,:,k)=-basis(2,:,k); 
     end
  end
    
  if Normal(2)==0
     if Normal(1)<0
        basis(2,:,k)=-basis(2,:,k); 
     end
  end
  
  if abs(Normal(1))==abs(Normal(2))
     if Normal(1)<0
        basis(2,:,k)=-basis(2,:,k); 
     end
  end

   Normal2=Normal;
   

    if Normal(1)==0
     if Normal(2)<0
       Normal2=-Normal;
       
     end
    end
    
  if Normal(2)==0
     if Normal(1)<0
       Normal2=-Normal;
       
     end
  end
  
  if abs(Normal(1))==abs(Normal(2))
     if Normal(1)<0
        Normal2=-Normal;
        
     end
  end 

normal(2,1,k)=Normal2(1);
normal(2,2,k)=Normal2(2); 
    
    %%%Third Edge Basis

 basis(3,2,k)=edge_length(new_ele(k,3))/(2*area);
 basis(3,1,k)=-basis(3,2,k)*X2;
 basis(3,3,k)=-basis(3,2,k)*Y2;
 
 
  
    midpoint=[(X3+X1)*0.5;(Y3+Y1)*0.5];
        
    Normal(1)=Y3-Y1;
    Normal(2)=X1-X3;
    Nlength=norm(Normal);
    Normal=Normal./Nlength;
   
    if distance2(midpoint+Normal,[X2;Y2])<distance2(midpoint-Normal,[X2;Y2])
        Normal=-Normal;
    end

  if Normal(1)==0
     if Normal(2)<0
        basis(3,:,k)=-basis(3,:,k); 
     end
  end
    
  if Normal(2)==0
     if Normal(1)<0
        basis(3,:,k)=-basis(3,:,k); 
     end
  end
  
  if abs(Normal(1))==abs(Normal(2))
     if Normal(1)<0
        basis(3,:,k)=-basis(3,:,k); 
     end
  end
  
   Normal2=Normal;
  
    %%Normal2 for prolongation
    if Normal(1)==0
     if Normal(2)<0
       Normal2=-Normal;
       
     end
    end
    
  if Normal(2)==0
     if Normal(1)<0
       Normal2=-Normal;
       
     end
  end
  
  if abs(Normal(1))==abs(Normal(2))
     if Normal(1)<0
        Normal2=-Normal;
        
     end
  end 

normal(3,1,k)=Normal2(1);
normal(3,2,k)=Normal2(2);
   


   
   basis(4,2,k)=1;

end