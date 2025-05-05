clear all; clc
%% Contour plots for Rg and EtE

L = 3.5*11; % length of fully open chain
N = 20; % number of residues
b = L/N; % Kuhn length

%% Storing coordinates of backbone beads
Backbone_data = xlsread('Backbone coordinates.csv'); % storing entire data

Data_size = size(Backbone_data); % size of data

No_Data_points = Data_size(1,1); % Number of data points

No_structures = No_Data_points/ N; % Number of frames

% storing individual structure

i = 1;
j = 1;

for p = 1:N: No_Data_points

         for j = 1:N
            Structure(i).x(j) = Backbone_data(p,1);
            Structure(i).y(j) = Backbone_data(p,2);
            Structure(i).z(j) = Backbone_data(p,3);
            p = p + 1; 
         end
         i = i+1; 

end


%% calculating Rg of each structure 

% calculating the COM coordinates of each structure

i = 1;
j = 1;

for  i = 1: No_structures
    xSum = 0;
    ySum = 0;
    zSum = 0;
    
    for j = 1:N
        xSum = xSum + Structure(i).x(j);
        ySum = ySum + Structure(i).y(j);
        zSum = zSum + Structure(i).z(j);
    end
    
    xCOM = xSum/N;
    yCOM = ySum/N;
    zCOM = zSum/N;
    
    j = 1;
     
    for j = 1:N
        xSquareDiff(j) = (Structure(i).x(j) - xCOM)^2;
        ySquareDiff(j) = (Structure(i).y(j) - yCOM)^2;
        zSquareDiff(j) = (Structure(i).z(j) - zCOM)^2;
    end
    
    Sum_Square_diff(i) = sum(xSquareDiff) + sum(ySquareDiff) + sum(zSquareDiff);
    avg_Sum_Square_diff(i) = Sum_Square_diff(i)/N;
    Rg(i)= sqrt(avg_Sum_Square_diff(i));
    
end


Rg = transpose(Rg);

%% storing coordinates of first and last beads of the backbone

% Backbone_ends_data = xlsread('Backbone end coordinates.xlsx'); % storing entire data
% 
% Data_size = size(Backbone_ends_data); % size of data
% 
% No_Data_points = Data_size(1,1); % Number of data points
% 
% i = 1;
% j = 1;
% 
% % storing coordinates of one end
% for i = 1:2:2*No_structures
%     End1x(j) = Backbone_ends_data(i,1);
%     End1y(j) = Backbone_ends_data(i,2);
%     End1z(j) = Backbone_ends_data(i,3);
%     j = j + 1;
% end
% 
% i = 1;
% j = 1; 
% 
% % storing coordinates of other end
% for i = 1:2:2*No_structures
%     End2x(j) = Backbone_ends_data(i+1,1);
%     End2y(j) = Backbone_ends_data(i+1,2);
%     End2z(j) = Backbone_ends_data(i+1,3);
%     j = j + 1;
% end
% 
% % calculating delx, dely, delz values
% 
% i = 1;
% for i = 1:No_structures
%     delx(i) = End2x(i)-End1x(i);
%     dely(i) = End2y(i)-End1y(i);
%     delz(i) = End2z(i)-End1z(i);
% end
% 
% delx = delx(:); 
% dely = dely(:);
% delz = delz(:);
% 
% EtE_dist = sqrt(delx.^2 + dely.^2 + delz.^2);

i = 1;
j = 1;


for i = 1:No_structures
    
    xDiff = (Structure(i).x(1)-Structure(i).x(N))^2;
    yDiff = (Structure(i).y(1)-Structure(i).y(N))^2;
    zDiff = (Structure(i).z(1)-Structure(i).z(N))^2;
    
    EtE_dist(i) = sqrt(xDiff + yDiff + zDiff);
    
end

EtE_dist = transpose(EtE_dist); 

%% plotting end to end distance
 j =1;
 time(1) = 0;
 
 for j = 1:No_structures-1
     time(j+1) = time(j) + 0.45237;
 end
 
 time = time(:);
 
figure(1)
plot(time,EtE_dist)
xlabel('Time (in ps)', 'FontSize', 40)
ylabel('End to End distance', 'FontSize', 40)
hold off



%% Define grid

min_EtE = min(EtE_dist);
max_EtE = max(EtE_dist);
min_Rg = min(Rg);
max_Rg = max(Rg);

% Define grid
xgrid = linspace(min_EtE, max_EtE);
ygrid = linspace(min_Rg, max_Rg);
[x1, y1] = meshgrid(xgrid,ygrid);

% Perfrom Kernel density estimate 
% xi is the desired grid points to evaluate
% f is an estimate of the density, eps give he EtE and Rg values
xi = [x1(:) y1(:)];
[ f, ep] = ksdensity([EtE_dist Rg], xi);


% format data in matrix for contourf and plot
X = reshape(ep(:,1),length(xgrid),length(ygrid));
Y = reshape(ep(:,2),length(xgrid),length(ygrid));
Z = reshape(f,length(xgrid),length(ygrid));
contourf(X,Y,Z,10)
ax=gca;
ax.XAxis.Exponent = 0;
xtickformat('%.4f')
xlabel('End to End Distance (in )')
ylabel('Radius of Gyration (in )')
colorbar

mean (Rg)
std(Rg)
mean(EtE_dist)
std(EtE_dist)
max(Z);
max(ans)