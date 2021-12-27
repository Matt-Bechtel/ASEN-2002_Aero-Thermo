close all
clear 
clc

format longG

T = 300;
P = 83500;
R = 287;

A2 = 1;
A1 = 9.5;

%[v_pitot v_venturi] = calcVelo();

    venturi_data = readmatrix('2002Data/Water_Manometer_data/venturi_water_data.csv');
    pitot_data = readmatrix('2002Data/Water_Manometer_data/pitot_water_data.csv');
    
    % sort data by voltage
    venturi_data = sort(venturi_data);
    pitot_data = sort(pitot_data);
    
    % convert to Pa
    venturi_data(:, 2) = venturi_data(:,2) * 284.84;  % in H2O to Pa
    pitot_data(:, 2) = pitot_data(:, 2) * 284.84; % in H20 to Pa
    
    vent_index = ones(20, 1);
    pitot_index = ones(20, 1);
 
    %indices corresponding to pressure differentials
    vent_index(2:end) = vent_index(2:end) .* find(diff(venturi_data(:,1)) > .2);
    pitot_index(2:end) = pitot_index(2:end) .* find(diff(pitot_data(:,1)) > .2);
    
    % preallocate
    v_venturi = zeros(length(vent_index), 2);
    v_pitot = zeros(length(pitot_index), 2);
   
    avg_venturi_data = zeros(length(vent_index), 2);
    avg_pitot_data = zeros(length(pitot_index), 2);
    
    %% Average data for venturi tube
    for i = 1:length(vent_index)
        
        if i == length(vent_index)
            avg_venturi_data(i, 1) = venturi_data(end, 1);
            avg_venturi_data(i, 2) = mean(venturi_data(vent_index(i) + 1:length(venturi_data), 2));
            
            break
        end
        
        avg_venturi_data(i, 1) = venturi_data(vent_index(i + 1));
        avg_venturi_data(i, 2) = mean(venturi_data(vent_index(i) + 1:vent_index(i + 1), 2));
        
        if i == 1
            avg_venturi_data(1, 2) = mean(venturi_data(1:vent_index(i + 1), 2));
        end
    end  
    
    
%% Average data for pitot tube
    for i = 1:length(pitot_index)
        
        if i == length(pitot_index)
            avg_pitot_data(i, 1) = pitot_data(end, 1);
            avg_pitot_data(i, 2) = mean(pitot_data(pitot_index(i) + 1:length(pitot_data), 2));
            
            break
        end
        
        avg_pitot_data(i, 1) = pitot_data(pitot_index(i + 1));
        avg_pitot_data(i, 2) = mean(pitot_data(pitot_index(i) + 1:pitot_index(i + 1), 2));
        
        if i == 1
            avg_pitot_data(1, 2) = mean(pitot_data(1:pitot_index(i + 1), 2));
        end
    end
    
    %% Velocity calculation
    for j = 1:length(vent_index)
        v_venturi(j, 1) = avg_venturi_data(j,1);
        v_venturi(j, 2) = sqrt((2 * avg_venturi_data(j, 2) * R * T) / (P * (1 - (A2 / A1) ^ 2))); % need A2 / A1, T, P
        
        v_pitot(j, 1) = avg_pitot_data(j,1);
        v_pitot(j, 2) = sqrt(2 * avg_pitot_data(j, 2) * ((R * T) / P)); % need T, P
    end
   
%% Plotting/misc things

    %avg_venturi_data(1,1) = 0.5;
    figure
    plot(venturi_data(:,1), venturi_data(:,2))
    hold on
    plot(pitot_data(:,1), pitot_data(:,2))
    xlabel('Voltage [V]')
    ylabel('Pressure Differential [Pa]')
    hold off
    legend('Venturi Data', 'Pitot Data')
    
    
    
    %% Velocity voltage plot
    figure
    plot(v_pitot(:, 1), v_pitot(:, 2))
    hold on
    plot(v_venturi(:, 1), v_venturi(:,2))
    xlabel('Voltage')
    ylabel('Velocity')
    legend('Pitot', 'Venturi')
    title('Velocity vs Voltage')
    
    
    
    
