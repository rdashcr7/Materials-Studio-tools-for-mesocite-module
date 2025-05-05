%% This code can be used to compare the structural deviation of a simulated mesomolecule from its actual structure's coordinates (mostly extracted from a PDB file) stored in an xlsx file.
%% Make sure to specify the number of monomer units 'N' of the mesomolecule (line 13). 
%% This code computes the radius of gyration of the simulated mesomolecule (referred to as backbone)
%% This code also computes Root Mean Square Deviation of the siimulated mesomolecule from the coordinates of actual structure (PDB structure.xlsx in this code)
%% This code also plots the probability contour plots of Rg and RMSD. You can also specify the number of desired contour lines in the contourf command.
%% The rotation matrix algorithm is used from the theory in 'https://faculty.sites.iastate.edu/jia/files/inline-files/rotation.pdf'.


clear all; clc

%% Radius of gyration and RMSD contour plots

N = 12; % number of residues

%% calculating Rg values from backbone coordinates

%% Storing coordinates of backbone beads
Backbone_data = xlsread('Coordinates.xlsx'); % storing entire data

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
%% storing coordinates of actual structure
Structure_data = xlsread('PDB Structure.xlsx'); % storing strcuture data

% shifting structure to origin
for i  = 1:N
    Structure_data_x(i) = Structure_data(i,1)-Structure_data(1,1);
    Structure_data_y(i) = Structure_data(i,2)-Structure_data(1,2);
    Structure_data_z(i) = Structure_data(i,3)-Structure_data(1,3);
end

Structure_data_x = Structure_data_x(:);
Structure_data_y = Structure_data_y(:);
Structure_data_z = Structure_data_z(:);

%% Storing coordinates of backbone beads
Backbone_data = xlsread('Coordinates.xlsx'); % storing entire data

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


%% shifting backbone trajectories to origin

i = 1;
j = 1;

for  i = 1 : No_structures
    basex(i) = Structure(i).x(1);
    basey(i) = Structure(i).y(1);
    basez(i) = Structure(i).z(1);
    for j = 1:N
        Structure(i).x(j) = Structure(i).x(j)-basex(i); % x cooordinate of first backbone bead is 0
        
        Structure(i).y(j) = Structure(i).y(j)-basey(i); % y cooordinate of first backbone bead is 0
        
        Structure(i).z(j) = Structure(i).z(j)-basez(i); % z cooordinate of first backbone bead is 0
    end
end

%% Rotation Matrix and rotaion

i = 1;
j = 1;

for  i = 1 : No_structures
    phi(i) = acosd((Structure(i).x(2)*Structure_data_x(2)+Structure(i).y(2)*Structure_data_y(2))/sqrt(((Structure(i).x(2))^2+(Structure(i).y(2))^2)*((Structure_data_x(2))^2+(Structure_data_y(2))^2))); %rotation angle phi
    theta(i) = acosd((Structure(i).x(2)*Structure_data_x(2)+Structure(i).z(2)*Structure_data_z(2))/sqrt(((Structure(i).x(2))^2+(Structure(i).z(2))^2)*((Structure_data_x(2))^2+(Structure_data_z(2))^2))); %rotation angle theta
    psi(i) = acosd((Structure(i).z(2)*Structure_data_z(2)+Structure(i).y(2)*Structure_data_y(2))/sqrt(((Structure(i).z(2))^2+(Structure(i).y(2))^2)*((Structure_data_z(2))^2+(Structure_data_y(2))^2))); %rotation angle psi
    Term_1(i) = cos(phi(i))*cos(theta(i));
    Term_2(i) = cos(phi(i))*sin(theta(i))*sin(psi(i))-sin(psi(i))*cos(psi(i));
    Term_3(i) = cos(phi(i))*sin(theta(i))*cos(psi(i))+sin(phi(i))*sin(psi(i));
    Term_4(i) = sin(phi(i))*cos(theta(i));
    Term_5(i) = sin(phi(i))*sin(theta(i))*sin(psi(i))+cos(phi(i))*cos(psi(i));
    Term_6(i) = sin(phi(i))*sin(theta(i))*cos(psi(i))-cos(phi(i))*sin(psi(i));
    Term_7(i) = -sin(theta(i));
    Term_8(i) = cos(theta(i))*sin(psi(i));
    Term_9(i) = cos(theta(i))*cos(psi(i));
    for j = 1:N
        Structure(i).x(j) = Structure(i).x(j)*Term_1(i) + Structure(i).y(j)*Term_2(i) + Structure(i).z(j)*Term_3(i);
        Structure(i).y(j) = Structure(i).x(j)*Term_4(i) + Structure(i).y(j)*Term_5(i) + Structure(i).z(j)*Term_6(i);
        Structure(i).z(j) = Structure(i).x(j)*Term_7(i) + Structure(i).y(j)*Term_8(i) + Structure(i).z(j)*Term_9(i);
    end
end

%% RMSD calculation

i = 1;

for i = 1: No_structures
    sum = 0;
    for j = 1:N
        xdiff(j) = (Structure(i).x(j)-Structure_data_x(j))^2;
        ydiff(j) = (Structure(i).y(j)-Structure_data_y(j))^2;
        zdiff(j) = (Structure(i).z(j)-Structure_data_z(j))^2;
        sum = sum + xdiff(j) + ydiff(j) +zdiff(j);
    end
    RMSD(i) = sqrt(sum/N);
end

RMSD = transpose(RMSD);
%% Define grid

min_RMSD = min(RMSD);
max_RMSD = max(RMSD);
min_Rg = min(Rg);
max_Rg = max(Rg);

% Define grid
xgrid = linspace(min_RMSD, max_RMSD);
ygrid = linspace(min_Rg, max_Rg);
[x1, y1] = meshgrid(xgrid,ygrid);

% Perfrom Kernel density estimate 
% xi is the desired grid points to evaluate
% f is an estimate of the density, eps give he EtE and Rg values
xi = [x1(:) y1(:)];
[ f, ep] = ksdensity([RMSD Rg], xi);

% format data in matrix for contourf and plot
X = reshape(ep(:,1),length(xgrid),length(ygrid));
Y = reshape(ep(:,2),length(xgrid),length(ygrid));
Z = reshape(f,length(xgrid),length(ygrid));
contourf(X,Y,Z,10)
ax=gca;
ax.XAxis.Exponent = 0;
xtickformat('%.4f')
xlabel('Root Mean Square Deviation (Å)','FontSize', 30,'FontWeight','bold')
ylabel('Radius of Gyration (Å)','FontSize', 30, 'FontWeight', 'bold')
title('Contour plot','FontSize', 30, 'FontWeight', 'bold')
colorbar
