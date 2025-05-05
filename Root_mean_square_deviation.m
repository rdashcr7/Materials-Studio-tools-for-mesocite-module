 clear all ; clc
 
 %% Calculating RMSD of structure
 
L = 3.5*11; % length of fully open chain
N = 12; % number of residues
b = L/N; % Kuhn length

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
Backbone_data = xlsread('Backbone coordinates.xlsx'); % storing entire data

Data_size = size(Backbone_data); % size of data

No_Data_points = Data_size(1,1); % Number of data points

No_structures = No_Data_points/ N; % Number of frames

% storing individual structure

i = 1;
j = 1;

for p = 1:12: No_Data_points

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

        
RMSD = RMSD.*0.1;

%% plotting RMSD vs time

%% plotting end to end distance
 j =1;
 time(1) = 0;
 
 for j = 1:No_Data_points/12-1
     time(j+1) = time(j) + 0.5104;
 end
 
 time = time(:);
 
figure(1)
plot(time,RMSD)
xlabel('Time (in ps)', 'FontSize', 40)
ylabel('RMSD (in nm)', 'FontSize', 40)
hold off