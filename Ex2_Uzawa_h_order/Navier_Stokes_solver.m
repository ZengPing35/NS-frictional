function [uh1,uh2,ph] = Navier_Stokes_solver(left,right,bottom,top,t_start,t_end,tau,h_partition,Gauss_point_number,basis_type_u,basis_type_p,max_iteration_step,tolerence,time_figure_err,figure_p,time_interval_p,figure_u,time_interval_u)
%h_partition: is the step size of the partition.
%basis_type: the type of the FE.
%basis_type_p=1:2D Lagrange linear FE.
%basis_type_p=11:2D P1b FE.
%basis_type_u=2:2D Lagrange quadratic FE.
%N1_basis,N2_basis:The N1 and N2 for the FE basis functions,not the partition.
%N1_partition,N2_partition:The N1 and N2 for the partition,not the FE basis functions.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.
%P_partition,T_partition, P_basis,T_basis: see the note in "generate_P_T_triangular,m".
%funciton_f1,funciton_f2: the right hand side functions of the Steady_Navier_Stokes equation.
%function_g1,function_g2: the Dirichelet boundary functions for (u1,u2)
%h_basis is the step size for the finite element nodes, not the partition.

%% Parameters
eps=10^(-8);   %The integral of p over the region is equal to 0
rho = 1;

%% The mesh
N1_partition = (right-left)/h_partition(1);
N2_partition = (top-bottom)/h_partition(2);
%  basis_type_u==11
u_number_basis = (N1_partition+1)*(N2_partition+1)+2*N1_partition*N2_partition;
uu_matrix_size = [u_number_basis u_number_basis];
u_vector_size = u_number_basis;
%  basis_type_p==1
p_N1_basis = N1_partition;
p_N2_basis = N2_partition;
p_number_basis = (p_N1_basis+1)*(p_N2_basis+1);
up_matrix_size=[u_number_basis p_number_basis];
pp_matrix_size = [p_number_basis p_number_basis];
p_vector_size=p_number_basis;
%  Mesh information for partition and finite element basis functions.
[P_partition,T_partition]=generate_P_T_triangle(left,right,bottom,top,h_partition,1);
%  basis_type_p==1
P_basis_p=P_partition;
T_basis_p=T_partition;
%  basis_type_u==11
[P_basis_u,T_basis_u]=generate_P_T_triangle(left,right,bottom,top,h_partition,11);

%% Guass quadrature's points and weights on the refenrece triangle and reference interval.
[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

%% Parameters for Matrix
number_of_elements=2*N1_partition*N2_partition;
%  basis_type_u==11
u_basis_index=[1 2 3 4];  
u_basis_index_length=length(u_basis_index);
%  basis_type_p==1
p_basis_index=[1 2 3];    
p_basis_index_length=length(p_basis_index);

%% Assemble the matrix.
% Bilinear terms (a_0 + b)
A1=assemble_matrix_2D('function_nu',P_partition,T_partition,T_basis_u,T_basis_u,number_of_elements,uu_matrix_size,u_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,u_basis_index_length,basis_type_u,1,0,u_basis_index_length,basis_type_u,1,0);
A2=assemble_matrix_2D('function_nu',P_partition,T_partition,T_basis_u,T_basis_u,number_of_elements,uu_matrix_size,u_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,u_basis_index_length,basis_type_u,0,1,u_basis_index_length,basis_type_u,0,1);
A3=assemble_matrix_2D('function_nu',P_partition,T_partition,T_basis_u,T_basis_u,number_of_elements,uu_matrix_size,u_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,u_basis_index_length,basis_type_u,1,0,u_basis_index_length,basis_type_u,0,1);
A5=assemble_matrix_2D('function_negativeone',P_partition,T_partition,T_basis_p,T_basis_u,number_of_elements,up_matrix_size,p_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,p_basis_index_length,basis_type_p,0,0,u_basis_index_length,basis_type_u,1,0);
A6=assemble_matrix_2D('function_negativeone',P_partition,T_partition,T_basis_p,T_basis_u,number_of_elements,up_matrix_size,p_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,p_basis_index_length,basis_type_p,0,0,u_basis_index_length,basis_type_u,0,1);
A7=eps*assemble_matrix_2D('function_negativeone',P_partition,T_partition,T_basis_p,T_basis_p,number_of_elements,pp_matrix_size,p_basis_index,p_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,p_basis_index_length,basis_type_p,0,0,p_basis_index_length,basis_type_p,0,0);
% O1=sparse(p_number_basis,p_number_basis);
A_bilinear=[2*A1+A2 A3 A5; A3' 2*A2+A1 A6;A5' A6' A7];
%  Assemble the mass matrix.(Term with time)
Me=assemble_matrix_2D('function_one',P_partition,T_partition,T_basis_u,T_basis_u,number_of_elements,uu_matrix_size,u_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,u_basis_index_length,basis_type_u,0,0,u_basis_index_length,basis_type_u,0,0);
O2=sparse(u_number_basis,u_number_basis);
O3=sparse(u_number_basis,p_number_basis);
O4=sparse(p_number_basis,p_number_basis);
M=[Me O2 O3;O2 Me O3;O3' O3' O4];

%% Assemble the load vector. 
% Source term
b1_f1=assemble_vector_2D_time('function_f1',N2_partition,P_partition,T_partition,T_basis_u,number_of_elements,u_vector_size,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type_u,0,0);
b2_f2=assemble_vector_2D_time('function_f2',N2_partition,P_partition,T_partition,T_basis_u,number_of_elements,u_vector_size,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type_u,0,0); 
O5=sparse(p_number_basis,1);
b_f=[b1_f1;b2_f2;O5];

%% Initialize the iteration in time.
u1_old_time=get_initial_vector('function_initial_u1',P_basis_u);
u2_old_time=get_initial_vector('function_initial_u2',P_basis_u);
p_old_time=get_initial_vector('function_initial_p',P_basis_p);
X_old_time=[u1_old_time;u2_old_time;p_old_time];
Lamda_old = zeros(2,(N1_partition-1));

%% Get the information matrices for boundary nodes and boundary edges.    
[boundary_nodes,boundary_edges]=generate_boundary_nodes_edges_Navier_Stokes(N1_partition,N2_partition,N1_partition,N2_partition);
[s,t] = find(boundary_nodes(1,:)==-2);
Lamda_old(1,:) = boundary_nodes(3,t);


%% Iteration in time.
N_t=(t_end-t_start)/tau;
for n=1:N_t            
    current_time=t_start+(tau*(n));    
%     fprintf('Time T = %d\n',current_time);

    %% Assemble the matrix of nonlinear part.
    AN2=assemble_matrix_2D_nonlinear('FE_solution_triangle_index',u1_old_time,u_basis_index,basis_type_u,0,0,P_partition,T_partition,T_basis_u,T_basis_u,number_of_elements,uu_matrix_size,u_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,u_basis_index_length,basis_type_u,1,0,u_basis_index_length,basis_type_u,0,0);
    AN3=assemble_matrix_2D_nonlinear('FE_solution_triangle_index',u2_old_time,u_basis_index,basis_type_u,0,0,P_partition,T_partition,T_basis_u,T_basis_u,number_of_elements,uu_matrix_size,u_basis_index,u_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,u_basis_index_length,basis_type_u,0,1,u_basis_index_length,basis_type_u,0,0);
    AN=[AN2+AN3 O2 O3;O2 AN2+AN3 O3;O3' O3' O4];
    A=M/tau+A_bilinear+AN;
    
    %% The load vector
    b_time = b_f + M/tau * X_old_time;

    %% Uzawa iteration
    % initial values
    u1h_old=X_old_time(1:u_number_basis);
    u2h_old=X_old_time(u_number_basis+1:2*u_number_basis);
    X_old_Uzawa=X_old_time;
    Lamda_old_Uzawa = Lamda_old;
    for l=1:max_iteration_step    
          
        
        %% Assemble the load vector. (The interface condition on \Gamma)
        b1_Gamma = zeros(u_number_basis,1);             
        b2_Gamma = assemble_vector_Gamma(h_partition, Lamda_old_Uzawa, current_time,u_vector_size, boundary_nodes, P_partition, T_partition, T_basis_u, Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle, u_basis_index, basis_type_u, 0, 0);
        b_Gamma = [b1_Gamma;b2_Gamma;O5]; 
        
         %% The load vector
        b = b_time + b_Gamma;
        
        %% Dealing with boundary conditions
        [A,b]=treat_Dirichlet_boundary_Navier_Stokes_time('function_g1','function_g2','function_g1_gamma',A,b,boundary_nodes,P_basis_u,u_number_basis);
%         [A,b]=fix_pressure_p_at_one_point('function_additional_condition_p',current_time,A,b,u_number_basis,p_number_basis,P_basis_p,NUM);
          
        %% the solution for(l-1)
        X= A\b;
        err = norm(X_old_Uzawa-X);
%         fprintf('Current time T = %d iteration %d error: %d\n',current_time,l,err);
        if h_partition(1) == figure_p
             %% Uzawa iteration error
            if current_time == time_figure_err
                Err_Uzawa(l) = err;
            end
        end
        
        if err < tolerence
            fprintf('Current time T = %d:iterate %d steps convergence!\n',current_time,l);
            break;
        end
        if l == max_iteration_step
            fprintf('Current time T = %d Not convergence!\n',current_time);
        end
        u1h_old=X(1:u_number_basis);
        u2h_old=X(u_number_basis+1:2*u_number_basis);
        X_old_Uzawa=X;
        
        %% Update Lamda
        Lamda_Uzawa = Lamda_old_Uzawa;
        Lamda_old_Uzawa = function_update_Lamda(P_partition, Lamda_Uzawa, u2h_old, rho, current_time, N1_partition, boundary_nodes);
    end   
    u1_old_time = u1h_old;
    u2_old_time = u2h_old;
    X_old_time=X_old_Uzawa;
    Lamda_old = Lamda_old_Uzawa;
    
    %% Figure
    %% Figure of u
    if h_partition(1) == figure_u
        if rem(current_time,time_interval_u)==0

            x=left:h_partition(1):right;
            y=bottom:h_partition(2):top;
            [x_new,y_new]=meshgrid(x,y);
            number_of_nodes_x=N1_partition+1;
            number_of_nodes_y=N2_partition+1;
            %The numerical solution of u
            u1=zeros(number_of_nodes_y,number_of_nodes_x);
            u2=zeros(number_of_nodes_y,number_of_nodes_x);
            for i=1:number_of_nodes_x
                u1(:,i)=u1h_old(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
                u2(:,i)=u2h_old(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
            end  

            figure;
            q = quiver(x_new,y_new,u1,u2);
            title('u');
            colorbar;
            colormap(jet);
            xlabel('x'),ylabel('y');
            mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), reshape(q.WData, numel(q.UData), [])).^2, 2));
            currentColormap = colormap(gca);
            [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
            cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            cmap(:,:,4) = 255;
            cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:3,:,:), [], 4).'); 
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:2,:,:), [], 4).');
            picturename = strcat('uh_',num2str(1/figure_u),'_time',num2str(current_time),'.fig');
            saveas(gca,picturename,'fig');
            close;
        end
    end
    %% Figure of p
    if h_partition(1) == figure_p
        
         %% Uzawa iteration error
        if current_time == time_figure_err
            filename = strcat('Err_Uzawa_time',num2str(current_time),'_iterate',num2str(l),'h_',num2str(1/h_partition(1)),'.mat');
            save(filename, 'Err_Uzawa');
            %% figure Uzawa iteration error
            figure;
%             plot(log10(Err_Uzawa(1,:)), '-.r','LineWidth',4);
            plot((Err_Uzawa(1,:)), '-.r^','LineWidth',4);
            xlabel('Iterations','interpreter','latex');
%             ylabel('$$\log(\|( u_1^{n, (m)},  u_2^{n, (m)}, p^{n, (m)}) - ( u_1^{n, (m-1)}, u_2^{n, (m-1)}, p^{n, (m-1)})\|_{L^2})$$','interpreter','latex');
            title('Error','interpreter','latex');
            a = 1/h_partition(1);
            picturename = strcat('Uzawa_iteration_error',num2str(a),'.fig');
            saveas(gca,picturename,'fig');
            close;
        end
        %% Figure
        if rem(current_time,time_interval_p)==0

            x=left:h_partition(1):right;
            y=bottom:h_partition(2):top;
            number_of_nodes_x=N1_partition+1;
            number_of_nodes_y=N2_partition+1;
            %The numerical solution of p
            p=zeros(number_of_nodes_y,number_of_nodes_x);
            ph_old=X_old_time(2*u_number_basis+1:2*u_number_basis+p_number_basis);
            for i=1:number_of_nodes_x
                p(:,i)=ph_old(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
            end  

            m = [0 1];
            qq = [-1 1];
            figure;
            imagesc(m,qq,p);
            shading interp;
            axis xy
            title('p');
            colorbar;
            colormap(jet);
            xlabel('x'),ylabel('y');
            picturename = strcat('ph_',num2str(1/figure_p),'_time',num2str(current_time),'.fig');
            saveas(gca,picturename,'fig');
            close;
        end
    end
end

%% Get the numerical solution for u1, u2 and p.
uh1=X(1:u_number_basis);
uh2=X(u_number_basis+1:2*u_number_basis);
ph=X(2*u_number_basis+1:2*u_number_basis+p_number_basis);