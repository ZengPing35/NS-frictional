function main

%% Parameters
format short e
%The problem domain is [left,right]*[bottom,top].
left=0;
right=1;
bottom=-1;
top=1;
basis_type_u=11;
basis_type_p=1;
Gauss_point_number=9;
max_iteration_step=5000;
tolerence=10^(-4);
t_start=0;
t_end=0.5;
figure_u = 1/16;
figure_p = 1/64;
time_interval_p = 0.05;
time_interval_u = 0.05;
time_figure_err = 0.05;

%% 保存输出
diary('data.txt')

%% Fix tau = 1/400
tau = 1/200

%% Exact solution （h = [1/64,1/64]）
h_partition=[1/64,1/64]                 
[uh1_exact,uh2_exact,ph_exact] = Navier_Stokes_solver(left,right,bottom,top,t_start,t_end,tau,h_partition,Gauss_point_number,basis_type_u,basis_type_p,max_iteration_step,tolerence,time_figure_err,figure_p,time_interval_p,figure_u,time_interval_u);
N2_partition_fine = (top-bottom)/h_partition(2);

filename = strcat('uh1_exact_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh1_exact');
filename = strcat('uh2_exact_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh2_exact');
filename = strcat('ph_exact_',num2str(1/h_partition(1)),'.mat');
save(filename, 'ph_exact');


%% Numerical solution (h = 1/4)
h_partition=[1/4,1/4]
[uh1,uh2,ph] = Navier_Stokes_solver(left,right,bottom,top,t_start,t_end,tau,h_partition,Gauss_point_number,basis_type_u,basis_type_p,max_iteration_step,tolerence,time_figure_err,figure_p,time_interval_p,figure_u,time_interval_u);

filename = strcat('uh1_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh1');
filename = strcat('uh2_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh2');
filename = strcat('ph_',num2str(1/h_partition(1)),'.mat');
save(filename, 'ph');

u_basis_index=[1 2 3 4]; 
p_basis_index=[1 2 3];

% L2-norm
uh1_L2_error=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
uh2_L2_error=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
ph_L2_error=FE_solution_error_triangle_index(ph,ph_exact,left,right,bottom,top,h_partition,p_basis_index,basis_type_p,0,0,Gauss_point_number,N2_partition_fine);
uh_L2_error_4=sqrt(uh1_L2_error^2+uh2_L2_error^2)
ph_L2_error_4=sqrt(ph_L2_error^2)

% H1-norm
uh1_H1_error_x=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh1_H1_error_y=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh1_H1_error=sqrt(uh1_H1_error_x^2+uh1_H1_error_y^2);
uh2_H1_error_x=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh2_H1_error_y=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh2_H1_error=sqrt(uh2_H1_error_x^2+uh2_H1_error_y^2);
uh_H1_error_4=sqrt(uh1_H1_error^2+uh2_H1_error^2)
uh_H1_4 = sqrt(uh_L2_error_4^2+uh_H1_error_4^2)

% Figure
x=left:h_partition(1):right;
y=bottom:h_partition(2):top;
N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);
[x_new,y_new]=meshgrid(x,y);
number_of_nodes_x=N1_partition+1;
number_of_nodes_y=N2_partition+1;
%The numerical solution of u and p
u1=zeros(number_of_nodes_y,number_of_nodes_x);
u2=zeros(number_of_nodes_y,number_of_nodes_x);
p=zeros(number_of_nodes_y,number_of_nodes_x);
for i=1:number_of_nodes_x
    u1(:,i)=uh1(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    u2(:,i)=uh2(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    p(:,i)=ph(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
end  

m = [0 1];
qq = [-1 1];
figure;
subplot(1,2,1);
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
subplot(1,2,2);
imagesc(m,qq,p);
shading interp;
axis xy
title('p');
colorbar;
colormap(jet);
xlabel('x'),ylabel('y');
picturename = strcat('u_p_t_end_4.fig');
saveas(gca,picturename,'fig');
close;


%% Numerical solution (h = 1/8)
h_partition=[1/8,1/8]
[uh1,uh2,ph] = Navier_Stokes_solver(left,right,bottom,top,t_start,t_end,tau,h_partition,Gauss_point_number,basis_type_u,basis_type_p,max_iteration_step,tolerence,time_figure_err,figure_p,time_interval_p,figure_u,time_interval_u);

filename = strcat('uh1_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh1');
filename = strcat('uh2_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh2');
filename = strcat('ph_',num2str(1/h_partition(1)),'.mat');
save(filename, 'ph');

u_basis_index=[1 2 3 4]; 
p_basis_index=[1 2 3];

% L2-norm
uh1_L2_error=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
uh2_L2_error=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
ph_L2_error=FE_solution_error_triangle_index(ph,ph_exact,left,right,bottom,top,h_partition,p_basis_index,basis_type_p,0,0,Gauss_point_number,N2_partition_fine);
uh_L2_error_8=sqrt(uh1_L2_error^2+uh2_L2_error^2)
r=(log(uh_L2_error_4)-log(uh_L2_error_8))/(log(2));
fprintf('Convergence order of uh_L2 = %f***\n',r);
ph_L2_error_8=sqrt(ph_L2_error^2)
r=(log(ph_L2_error_4)-log(ph_L2_error_8))/(log(2));
fprintf('Convergence order of ph_L2 = %f***\n',r);

% H1-norm
uh1_H1_error_x=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh1_H1_error_y=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh1_H1_error=sqrt(uh1_H1_error_x^2+uh1_H1_error_y^2);
uh2_H1_error_x=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh2_H1_error_y=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh2_H1_error=sqrt(uh2_H1_error_x^2+uh2_H1_error_y^2);
uh_H1_error_8=sqrt(uh1_H1_error^2+uh2_H1_error^2)
r=(log(uh_H1_error_4)-log(uh_H1_error_8))/(log(2));
fprintf('Convergence order of uh_sime_H1 = %f***\n',r)
uh_H1_8 = sqrt(uh_L2_error_8^2+uh_H1_error_8^2)
r=(log(uh_H1_4)-log(uh_H1_8))/(log(2));
fprintf('Convergence order of uh_H1 = %f***\n',r)

% Figure
x=left:h_partition(1):right;
y=bottom:h_partition(2):top;
N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);
[x_new,y_new]=meshgrid(x,y);
number_of_nodes_x=N1_partition+1;
number_of_nodes_y=N2_partition+1;
%The numerical solution of u and p
u1=zeros(number_of_nodes_y,number_of_nodes_x);
u2=zeros(number_of_nodes_y,number_of_nodes_x);
p=zeros(number_of_nodes_y,number_of_nodes_x);
for i=1:number_of_nodes_x
    u1(:,i)=uh1(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    u2(:,i)=uh2(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    p(:,i)=ph(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
end  
m = [0 1];
qq = [-1 1];
figure;
subplot(1,2,1);
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
subplot(1,2,2);
imagesc(m,qq,p);
shading interp;
axis xy
title('p');
colorbar;
colormap(jet);
xlabel('x'),ylabel('y');
picturename = strcat('u_p_t_end_8.fig');
saveas(gca,picturename,'fig');
close;


%% Numerical solution (h = 1/16)
h_partition=[1/16,1/16]
[uh1,uh2,ph] = Navier_Stokes_solver(left,right,bottom,top,t_start,t_end,tau,h_partition,Gauss_point_number,basis_type_u,basis_type_p,max_iteration_step,tolerence,time_figure_err,figure_p,time_interval_p,figure_u,time_interval_u);

filename = strcat('uh1_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh1');
filename = strcat('uh2_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh2');
filename = strcat('ph_',num2str(1/h_partition(1)),'.mat');
save(filename, 'ph');

u_basis_index=[1 2 3 4]; 
p_basis_index=[1 2 3];

% L2-norm
uh1_L2_error=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
uh2_L2_error=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
ph_L2_error=FE_solution_error_triangle_index(ph,ph_exact,left,right,bottom,top,h_partition,p_basis_index,basis_type_p,0,0,Gauss_point_number,N2_partition_fine);
uh_L2_error_16=sqrt(uh1_L2_error^2+uh2_L2_error^2)
r=(log(uh_L2_error_8)-log(uh_L2_error_16))/(log(2));
fprintf('Convergence order of uh_L2 = %f***\n',r);
ph_L2_error_16=sqrt(ph_L2_error^2)
r=(log(ph_L2_error_8)-log(ph_L2_error_16))/(log(2));
fprintf('Convergence order of ph_L2 = %f***\n',r);

% H1-norm
uh1_H1_error_x=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh1_H1_error_y=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh1_H1_error=sqrt(uh1_H1_error_x^2+uh1_H1_error_y^2);
uh2_H1_error_x=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh2_H1_error_y=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh2_H1_error=sqrt(uh2_H1_error_x^2+uh2_H1_error_y^2);
uh_H1_error_16=sqrt(uh1_H1_error^2+uh2_H1_error^2)
r=(log(uh_H1_error_8)-log(uh_H1_error_16))/(log(2));
fprintf('Convergence order of uh_sime_H1 = %f***\n',r)
uh_H1_16 = sqrt(uh_L2_error_16^2+uh_H1_error_16^2)
r=(log(uh_H1_8)-log(uh_H1_16))/(log(2));
fprintf('Convergence order of uh_H1 = %f***\n',r)

% Figure
x=left:h_partition(1):right;
y=bottom:h_partition(2):top;
N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);
[x_new,y_new]=meshgrid(x,y);
number_of_nodes_x=N1_partition+1;
number_of_nodes_y=N2_partition+1;
%The numerical solution of u and p
u1=zeros(number_of_nodes_y,number_of_nodes_x);
u2=zeros(number_of_nodes_y,number_of_nodes_x);
p=zeros(number_of_nodes_y,number_of_nodes_x);
for i=1:number_of_nodes_x
    u1(:,i)=uh1(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    u2(:,i)=uh2(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    p(:,i)=ph(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
end  
m = [0 1];
qq = [-1 1];
figure;
subplot(1,2,1);
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
subplot(1,2,2);
imagesc(m,qq,p);
shading interp;
axis xy
title('p');
colorbar;
colormap(jet);
xlabel('x'),ylabel('y');
picturename = strcat('u_p_t_end_16.fig');
saveas(gca,picturename,'fig');
close;


%% Numerical solution (h = 1/32)
h_partition=[1/32,1/32]
[uh1,uh2,ph] = Navier_Stokes_solver(left,right,bottom,top,t_start,t_end,tau,h_partition,Gauss_point_number,basis_type_u,basis_type_p,max_iteration_step,tolerence,time_figure_err,figure_p,time_interval_p,figure_u,time_interval_u);

filename = strcat('uh1_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh1');
filename = strcat('uh2_',num2str(1/h_partition(1)),'.mat');
save(filename, 'uh2');
filename = strcat('ph_',num2str(1/h_partition(1)),'.mat');
save(filename, 'ph');

u_basis_index=[1 2 3 4]; 
p_basis_index=[1 2 3];

% L2-norm
uh1_L2_error=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
uh2_L2_error=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
ph_L2_error=FE_solution_error_triangle_index(ph,ph_exact,left,right,bottom,top,h_partition,p_basis_index,basis_type_p,0,0,Gauss_point_number,N2_partition_fine);
uh_L2_error_32=sqrt(uh1_L2_error^2+uh2_L2_error^2)
r=(log(uh_L2_error_16)-log(uh_L2_error_32))/(log(2));
fprintf('Convergence order of uh_L2 = %f***\n',r);
ph_L2_error_32=sqrt(ph_L2_error^2)
r=(log(ph_L2_error_16)-log(ph_L2_error_32))/(log(2));
fprintf('Convergence order of ph_L2 = %f***\n',r);

% H1-norm
uh1_H1_error_x=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh1_H1_error_y=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh1_H1_error=sqrt(uh1_H1_error_x^2+uh1_H1_error_y^2);
uh2_H1_error_x=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
uh2_H1_error_y=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
uh2_H1_error=sqrt(uh2_H1_error_x^2+uh2_H1_error_y^2);
uh_H1_error_32=sqrt(uh1_H1_error^2+uh2_H1_error^2)
r=(log(uh_H1_error_16)-log(uh_H1_error_32))/(log(2));
fprintf('Convergence order of uh_sime_H1 = %f***\n',r)
uh_H1_32 = sqrt(uh_L2_error_32^2+uh_H1_error_32^2)
r=(log(uh_H1_16)-log(uh_H1_32))/(log(2));
fprintf('Convergence order of uh_H1 = %f***\n',r)

% Figure
x=left:h_partition(1):right;
y=bottom:h_partition(2):top;
N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);
[x_new,y_new]=meshgrid(x,y);
number_of_nodes_x=N1_partition+1;
number_of_nodes_y=N2_partition+1;
%The numerical solution of u and p
u1=zeros(number_of_nodes_y,number_of_nodes_x);
u2=zeros(number_of_nodes_y,number_of_nodes_x);
p=zeros(number_of_nodes_y,number_of_nodes_x);
for i=1:number_of_nodes_x
    u1(:,i)=uh1(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    u2(:,i)=uh2(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
    p(:,i)=ph(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
end  
m = [0 1];
qq = [-1 1];
figure;
subplot(1,2,1);
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
subplot(1,2,2);
imagesc(m,qq,p);
shading interp;
axis xy
title('p');
colorbar;
colormap(jet);
xlabel('x'),ylabel('y');
picturename = strcat('u_p_t_end_32.fig');
saveas(gca,picturename,'fig');
close;


% %% Numerical solution (h = 1/64)
% h_partition=[1/64,1/64]
% [uh1,uh2,ph] = Navier_Stokes_solver(left,right,bottom,top,t_start,t_end,tau,h_partition,Gauss_point_number,basis_type_u,basis_type_p,max_iteration_step,tolerence,time_figure_err,figure_p,time_interval_p,figure_u,time_interval_u);
% 
% filename = strcat('uh1_',num2str(1/h_partition(1)),'.mat');
% save(filename, 'uh1');
% filename = strcat('uh2_',num2str(1/h_partition(1)),'.mat');
% save(filename, 'uh2');
% filename = strcat('ph_',num2str(1/h_partition(1)),'.mat');
% save(filename, 'ph');
% 
% u_basis_index=[1 2 3 4]; 
% p_basis_index=[1 2 3];
% 
% % L2-norm
% uh1_L2_error=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
% uh2_L2_error=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,0,Gauss_point_number,N2_partition_fine);
% ph_L2_error=FE_solution_error_triangle_index(ph,ph_exact,left,right,bottom,top,h_partition,p_basis_index,basis_type_p,0,0,Gauss_point_number,N2_partition_fine);
% uh_L2_error_64=sqrt(uh1_L2_error^2+uh2_L2_error^2)
% r=(log(uh_L2_error_32)-log(uh_L2_error_64))/(log(2));
% fprintf('Convergence order of uh_L2 = %f***\n',r);
% ph_L2_error_64=sqrt(ph_L2_error^2)
% r=(log(ph_L2_error_32)-log(ph_L2_error_64))/(log(2));
% fprintf('Convergence order of ph_L2 = %f***\n',r);
% 
% % H1-norm
% uh1_H1_error_x=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
% uh1_H1_error_y=FE_solution_error_triangle_index(uh1,uh1_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
% uh1_H1_error=sqrt(uh1_H1_error_x^2+uh1_H1_error_y^2);
% uh2_H1_error_x=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,1,0,Gauss_point_number,N2_partition_fine);
% uh2_H1_error_y=FE_solution_error_triangle_index(uh2,uh2_exact,left,right,bottom,top,h_partition,u_basis_index,basis_type_u,0,1,Gauss_point_number,N2_partition_fine);
% uh2_H1_error=sqrt(uh2_H1_error_x^2+uh2_H1_error_y^2);
% uh_H1_error_64=sqrt(uh1_H1_error^2+uh2_H1_error^2)
% r=(log(uh_H1_error_32)-log(uh_H1_error_64))/(log(2));
% fprintf('Convergence order of uh_sime_H1 = %f***\n',r)
% uh_H1_64 = sqrt(uh_L2_error_64^2+uh_H1_error_64^2)
% r=(log(uh_H1_32)-log(uh_H1_64))/(log(2));
% fprintf('Convergence order of uh_H1 = %f***\n',r)
% 
% % Figure
% x=left:h_partition(1):right;
% y=bottom:h_partition(2):top;
% N1_partition=(right-left)/h_partition(1);
% N2_partition=(top-bottom)/h_partition(2);
% [x_new,y_new]=meshgrid(x,y);
% number_of_nodes_x=N1_partition+1;
% number_of_nodes_y=N2_partition+1;
% %The numerical solution of u and p
% u1=zeros(number_of_nodes_y,number_of_nodes_x);
% u2=zeros(number_of_nodes_y,number_of_nodes_x);
% p=zeros(number_of_nodes_y,number_of_nodes_x);
% for i=1:number_of_nodes_x
%     u1(:,i)=uh1(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
%     u2(:,i)=uh2(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
%     p(:,i)=ph(number_of_nodes_y*(i-1)+1:number_of_nodes_y*i);
% end  
% m = [0 1];
% qq = [-1 1];
% figure;
% subplot(1,2,1);
% q = quiver(x_new,y_new,u1,u2);
% title('u');
% colorbar;
% colormap(jet);
% xlabel('x'),ylabel('y');
% mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), reshape(q.WData, numel(q.UData), [])).^2, 2));
% currentColormap = colormap(gca);
% [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
% cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
% cmap(:,:,4) = 255;
% cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
% set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:3,:,:), [], 4).'); 
% set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:2,:,:), [], 4).');
% subplot(1,2,2);
% imagesc(m,qq,p);
% shading interp;
% axis xy
% title('p');
% colorbar;
% colormap(jet);
% xlabel('x'),ylabel('y');
% picturename = strcat('u_p_t_end_64.fig');
% saveas(gca,picturename,'fig');
% close;

%% 保存输出
diary off;