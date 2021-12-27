clear all
close all
clc
format longG

%% directory load for pitot tube
pitotfiles = dir(['Aero Lab Windtunnel Calibration/Aero Lab 1 - 2019 Group Data/VelocityVoltageData/PitotProbeToPressureTransducer/*.csv']);

for i = 1:numel(pitotfiles);
    Pitot = load(['Aero Lab Windtunnel Calibration/Aero Lab 1 - 2019 Group Data/VelocityVoltageData/PitotProbeToPressureTransducer/',pitotfiles(i).name]);
    y = Pitot(:, 7);
    y_ind = find(diff(y)>1);
    y_ind = [0; y_ind; length(y)];
    
    %parses data based off index change from voltage to voltage
    PitotData = NaN(length(y_ind)-1, 7);
    for b = 1:(length(y_ind)-1);
        temp = mean(Pitot(y_ind(b) + 1:y_ind(b+1),:));
        PitotData(b,:) = temp;
    end
    
    %% Create struct to store data
    VenturiTrialHolder(i).v = PitotData;
    
for i = 1:4
    fivebyseven = NaN(5,7);
    
    for j = 1:5
        temporary = mean([VenturiTrialHolder(i).v(j,:);VenturiTrialHolder(i+4).v(j,:);VenturiTrialHolder(i+8).v(j,:)]);
        fivebyseven(j,:) = temporary;
    end
    
    CompressedVenturiTrialHolder(i).v = fivebyseven;
end

for i = 1:4
    fivebyseven = NaN(5,7);
    for j = 1:5
        temporary = mean([PitotTrialHolder(i).v(j,:);PitotTrialHolder(i+4).v(j,:);PitotTrialHolder(i+8).v(j,:)]);
        fivebyseven(j,:) = temporary;
    end
    CompressedPitotTrialHolder(i).v = fivebyseven;
end
    
end

clear i j
%% directory load for venturi tube
venturifiles = dir(['Aero Lab Windtunnel Calibration/Aero Lab 1 - 2019 Group Data/VelocityVoltageData/VenturiTubeToPressureTransducer/*.csv']);

for i = 1:numel(venturifiles);
load(['Aero Lab Windtunnel Calibration/Aero Lab 1 - 2019 Group Data/VelocityVoltageData/VenturiTubeToPressureTransducer/', venturifiles(i).name]);
end

%% parse out data

