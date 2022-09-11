%% Setup
clear
close all
clc

%% Go through the folders
folderList = dir('Data/');
folderList = {folderList(4:end).name};

outputTable = NaN.*ones(35, 2, 3); % Participant x Period x Movement

for folderNum = 1:length(folderList)
    folderName = folderList{folderNum}
    
    filesList = dir(strcat('Data/', folderName));
    filesList = {filesList(end-1:end).name};
    
    fileAData = readtable(strcat('Data/', folderName, '/', filesList{1}));
    fileAData = fileAData(~any(ismissing(fileAData), 2), :);
    fileBData = readtable(strcat('Data/', folderName, '/', filesList{2}));
    fileBData = fileBData(~any(ismissing(fileBData), 2), :);
    
    fileAData.Properties.VariableNames = {'Dev', 'EF', 'Rot'};
    fileBData.Properties.VariableNames = {'Dev', 'EF', 'Rot'};
    
    participantNum = strsplit(folderName, 'BID');
    participantNum = str2double(participantNum{2});
    
    periodA = strsplit(filesList{1}, 'Per');
    periodA = strsplit(periodA{2}, '.txt');
    periodA = str2double(periodA{1});
    
    periodB = strsplit(filesList{2}, 'Per');
    periodB = strsplit(periodB{2}, '.txt');
    periodB = str2double(periodB{1});
    
    outputTable(participantNum, periodA, 1) = ratioOfActiveTime(fileAData.Dev);
    outputTable(participantNum, periodA, 2) = ratioOfActiveTime(fileAData.EF);
    outputTable(participantNum, periodA, 3) = ratioOfActiveTime(fileAData.Rot);
    
    outputTable(participantNum, periodB, 1) = ratioOfActiveTime(fileBData.Dev);
    outputTable(participantNum, periodB, 2) = ratioOfActiveTime(fileBData.EF);
    outputTable(participantNum, periodB, 3) = ratioOfActiveTime(fileBData.Rot);
end

DataTable = mulDimArray2table(outputTable, {'Participant', 'Period', 'Movement'}, 'Value');