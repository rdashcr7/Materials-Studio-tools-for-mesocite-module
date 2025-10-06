clear all; clc

%% normalized 1-D gaussian PDF in each dimension 

N = input('Enter the number of amino acids in the peptide:'); % number of residues
L = 3.5*11; % length of fully open chain
b = L/N; % Kuhn length

% storing coordinates of first and last beads of the backbone

Backbone_ends_data = xlsread('Backbone Ends.xlsx'); % storing entire data

Data_size = size(Backbone_ends_data); % size of data

No_Data_points = Data_size(1,1); % Number of data points

i = 1;
j = 1;

% storing coordinates of one end
for i = 1:2:No_Data_points
    End1x(j) = Backbone_ends_data(i,1);
    End1y(j) = Backbone_ends_data(i,2);
    End1z(j) = Backbone_ends_data(i,3);
    j = j + 1;
end

i = 1;
j = 1; 

% storing coordinates of other end
for i = 1:2:No_Data_points
    End2x(j) = Backbone_ends_data(i+1,1);
    End2y(j) = Backbone_ends_data(i+1,2);
    End2z(j) = Backbone_ends_data(i+1,3);
    j = j + 1;
end

% calculating delx, dely, delz values

i = 1;
for i = 1:No_Data_points/2
    delx(i) = End2x(i)-End1x(i);
    dely(i) = End2y(i)-End1y(i);
    delz(i) = End2z(i)-End1z(i);
end

delx = delx(:); 
dely = dely(:);
delz = delz(:);

%% calculating 1-D probability distribution in x direction
avg_x_sq = mean(delx.*delx);
bin_width = 0.025; % in angstrom

Total_data = size(delx);
max_val = max(delx); % max end to end distance value
min_val = min(delx); % min end to end distance value 
delx_iter_val = min_val; % EtE value to be checked 
No_of_bins = round((max_val - min_val)/bin_width); % number of bins

bin =zeros(No_of_bins,1);

for j = 1 : No_of_bins
    for i = 1:Total_data
        if delx(i) >= delx_iter_val && delx(i) <= delx_iter_val + bin_width
            bin(j)= bin(j) + 1 ; 
        end
    end
     
    delx_iter_val = delx_iter_val + bin_width;
end

% calculating probability from simulation data
delx_initial = min_val;
for  j = 1:No_of_bins
    delx_initial = delx_initial + bin_width/2;
    delx_x_axis(j) = delx_initial;
    Probabilityx(j) = bin(j)/sum(bin)*(avg_x_sq^0.5)/bin_width;
%     Probabilityx(j) = bin(j)/sum(bin)/bin_width;
    delx_initial = delx_initial + bin_width/2;
end

% calculating probability from the formula
i = 1;
for i = 1:No_Data_points/2
    Px(i) = sqrt(3)/sqrt(2*pi*N*b*b)*exp(-3*(delx(i))^2/(2*N*b*b))*sqrt(avg_x_sq);
%     Px(i) = 1/sqrt(2*pi*avg_x_sq)*exp(-(delx(i))^2/(2*avg_x_sq));
end

% calculating probability by fitting data into normal distribution
Px_dist = fitdist(delx, 'Normal');
Px_pd = pdf(Px_dist, delx)*(avg_x_sq)^(0.5);
% Px_pd = pdf(Px_dist, delx);

% plotting
figure(1)
plot(delx_x_axis/avg_x_sq^0.5, Probabilityx, 'LineWidth', 1)
% plot(delx_x_axis, Probabilityx,'.','MarkerSize',1)
hold on
plot(delx/avg_x_sq^0.5, Px,'.', 'MarkerSize',10)
% plot(delx, Px,'.', 'MarkerSize',10)
hold on
plot(delx/avg_x_sq^0.5, Px_pd,'o')
legend({'Probability calculated from DPD', 'Probability calculated from function','Fitting data into a normal distribution'},'FontSize', 20)
% legend({'Probability calculated from DPD', 'Probability calculated from function'},'FontSize', 20)
xlabel('x/<x^{2}>^{1/2}', 'FontSize', 40)
ylabel('P_{1d}(N,x)<x^2>^{1/2}', 'FontSize', 40)
hold off
%% calculating 1-D probability distribution in y direction
avg_y_sq = mean(dely.*dely);

Total_data = size(dely);
max_val = max(dely); % max end to end distance value
min_val = min(dely); % min end to end distance value 
dely_iter_val = min_val; % EtE value to be checked 
No_of_bins = round((max_val - min_val)/bin_width); % number of bins

bin =zeros(No_of_bins,1);

for j = 1 : No_of_bins
    for i = 1:Total_data
        if dely(i) >= dely_iter_val && dely(i) <= dely_iter_val + bin_width
            bin(j)= bin(j) + 1 ; 
        end
    end
     
    dely_iter_val = dely_iter_val + bin_width;
end

% calculating probability from simulation data
dely_initial = min_val;
for  j = 1:No_of_bins
    dely_initial = dely_initial + bin_width/2;
    dely_x_axis(j) = dely_initial;
    Probabilityy(j) = bin(j)/sum(bin)*(avg_y_sq^0.5)/bin_width;
%     Probability(j) = bin(j)/sum(bin);
    dely_initial = dely_initial + bin_width/2;
end

% calculating probability from the formula
i = 1;
for i = 1:No_Data_points/2
    Py(i) = sqrt(3)/sqrt(2*pi*N*b*b)*exp(-3*(dely(i))^2/(2*N*b*b))*sqrt(avg_y_sq);
    %Py(i) = 1/sqrt(2*pi*avg_y_sq)*exp(-(dely(i))^2/(2*avg_y_sq));
end

% calculating probability by fitting data into normal distribution
Py_dist = fitdist(dely, 'Normal');
Py_pd = pdf(Py_dist, dely)*(avg_y_sq)^(0.5);

figure(2)
plot(dely_x_axis/avg_y_sq^0.5, Probabilityy, 'LineWidth', 1)
hold on
plot(dely/avg_y_sq^0.5, Py,'.', 'MarkerSize',10)
hold on
plot(dely/avg_y_sq^0.5, Py_pd,'o')
legend({'Probability calculated from DPD', 'Probability calculated from function','Fitting data into a normal distribution'},'FontSize', 20)
xlabel('y/<y^{2}>^{1/2}', 'FontSize', 40)
ylabel('P_{1d}(N,y)<y^2>^{1/2}', 'FontSize', 40)
hold off

%% calculating 1-D probability distribution in z direction
avg_z_sq = mean(delz.*delz);

Total_data = size(delz);
max_val = max(delz); % max end to end distance value
min_val = min(delz); % min end to end distance value 
delz_iter_val = min_val; % EtE value to be checked 
No_of_bins = round((max_val - min_val)/bin_width); % number of bins

bin =zeros(No_of_bins,1);

for j = 1 : No_of_bins
    for i = 1:Total_data
        if delz(i) >= delz_iter_val && delz(i) <= delz_iter_val + bin_width
            bin(j)= bin(j) + 1 ; 
        end
    end
     
    delz_iter_val = delz_iter_val + bin_width;
end

% calculating probability from simulation data
delz_initial = min_val;
for  j = 1:No_of_bins
    delz_initial = delz_initial + bin_width/2;
    delz_x_axis(j) = delz_initial;
    Probabilityz(j) = bin(j)/sum(bin)*(avg_z_sq^0.5)/bin_width;
%     Probabilityz(j) = bin(j)/sum(bin);
    delz_initial = delz_initial + bin_width/2;
end

% calculating probability from the formula
i = 1;
for i = 1:No_Data_points/2
    Pz(i) = sqrt(3)/sqrt(2*pi*N*b*b)*exp(-3*(delz(i))^2/(2*N*b*b))*sqrt(avg_z_sq);
    %Pz(i) = 1/sqrt(2*pi*avg_z_sq)*exp(-(delz(i))^2/(2*avg_z_sq));
end

% calculating probability by fitting data into normal distribution
Pz_dist = fitdist(delz, 'Normal');
Pz_pd = pdf(Pz_dist, delz)*(avg_z_sq)^(0.5);

figure(3)
plot(delz_x_axis/avg_z_sq^0.5, Probabilityz, 'LineWidth', 1)
hold on
plot(delz/avg_z_sq^0.5, Pz,'.', 'MarkerSize',10)
hold on
plot(delz/avg_z_sq^0.5, Pz_pd,'o')
legend({'Probability calculated from DPD', 'Probability calculated from function','Fitting data into a normal distribution'},'FontSize', 20)
xlabel('z/<z^{2}>^{1/2}', 'FontSize', 40)
ylabel('P_{1d}(N,z)<z^2>^{1/2}', 'FontSize', 40)
hold off
