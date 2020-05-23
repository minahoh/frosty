%%%%%%%%%%%%%%%%%%%%%
%Author: Minah Oh
%Date: 04/27/2020
%%%%%%%%%%%%%%%%%%%%%

function v=mg(level, u0, b, n, Qcell, Acell,VPatchCell1, VPatchCell2, New_ele, Edge) %u0 is the initial value

A=cell2mat(Acell(level));
Q=cell2mat(Qcell(level));

if level==1
   v=A\b;
end


if level>1

   
    Vpatch1=cell2mat(VPatchCell1(level));
    Vpatch2=cell2mat(VPatchCell2(level));

    new_ele=cell2mat(New_ele(level-1));
    edge=cell2mat(Edge(level-1));
    [N_ele,~]=size(new_ele);
    [N_edge,~]=size(edge);
    
    v=gs(A,b,u0,Vpatch1,Vpatch2);
  
    %Coarse grid correction
    zero_vec=zeros(N_edge+N_ele,1);
    diff=b-A*v;
    OP2=Q'*diff;
    v_short=mg(level-1, zero_vec,OP2,n,Qcell,Acell, VPatchCell1, VPatchCell2, New_ele, Edge);
    v=v+Q*v_short;
   
    v=gs_transpose(A,b,v,Vpatch1,Vpatch2);
  
end