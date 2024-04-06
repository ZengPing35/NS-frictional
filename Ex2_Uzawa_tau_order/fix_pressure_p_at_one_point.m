function [A,b]=fix_pressure_p_at_one_point(Dirichlet_boundary_function_name_p,current_time,A,b,number_of_FE_nodes_u,number_of_FE_nodes_p,P_basis_p,NUM)
% p without any derivative :addition condition for p
% NUM<number_of_FE_nodes_p
for k=1:number_of_FE_nodes_p    
    if k==NUM        
       p_index=2*number_of_FE_nodes_u+NUM;
       A(p_index,:)=0;
       A(p_index,p_index)=1;
       b(p_index,1)=feval(Dirichlet_boundary_function_name_p,P_basis_p(1,NUM),P_basis_p(2,NUM),current_time);
    end    
end


