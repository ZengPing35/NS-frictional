function figure_err

load('Err_Uzawa_time0.05_iterate43h_64.mat');
figure;
plot(log10(Err_Uzawa(1,:)), '-r^','LineWidth',4);
xlabel('Iterations','interpreter','latex');
ylabel('$$\log(E^m)$$','interpreter','latex');
title('Error','interpreter','latex');
picturename = strcat('Uzawa_iteration_error_64.fig');
saveas(gca,picturename,'fig');
close;