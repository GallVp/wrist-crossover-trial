%% Setup
clear
close all
clc

%% Parameters
FS = 1000;
fs = 25;

%% Load Data
HM_AE       = load('Data Files/HM AE');
HM_PA       = load('Data Files/HM PA');
HM_Seated   = load('Data Files/NWM_seated');

%% Low-pass filter
HM_Seated_roll  = lowPassStream(HM_Seated.data.roll,   FS, 10);
HM_Seated_pitch = lowPassStream(HM_Seated.data.pitch,  FS, 10);
HM_Seated_yaw   = lowPassStream(HM_Seated.data.yaw,    FS, 10);

%% Downsample
HM_Seated_roll_ds   = resample(HM_Seated_roll, 1, 40);
HM_Seated_pitch_ds  = resample(HM_Seated_pitch, 1, 40);
HM_Seated_yaw_ds    = resample(HM_Seated_yaw, 1, 40);


%% Find velocity
V_HM_Seated_roll    = diff(HM_Seated_roll_ds)./(1/fs);
V_HM_Seated_pitch   = diff(HM_Seated_pitch_ds)./(1/fs);
V_HM_Seated_yaw     = diff(HM_Seated_yaw_ds)./(1/fs);

%% Remove first and last second of data
V_HM_Seated_roll    = V_HM_Seated_roll(fs+1:end-fs);
V_HM_Seated_pitch   = V_HM_Seated_pitch(fs+1:end-fs);
V_HM_Seated_yaw     = V_HM_Seated_yaw(fs+1:end-fs);

%% Find thresholds
std(V_HM_Seated_roll)
std(V_HM_Seated_pitch)
std(V_HM_Seated_yaw)

%% Plot Data
ratioOfActiveTime(HM_AE.data.roll, 'HM AE roll');
ratioOfActiveTime(HM_AE.data.pitch, 'HM AE pitch');
ratioOfActiveTime(HM_AE.data.yaw, 'HM AE yaw');

ratioOfActiveTime(HM_PA.data.roll, 'HM PA roll');
ratioOfActiveTime(HM_PA.data.pitch, 'HM PA pitch');
ratioOfActiveTime(HM_PA.data.yaw, 'HM PA yaw');