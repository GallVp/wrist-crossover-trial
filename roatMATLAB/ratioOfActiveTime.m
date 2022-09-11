function r = ratioOfActiveTime(signal, titleText)

if nargin < 2
    titleText = '';
end

%% Parameters
FS                  = 1000;
FS_DS               = 25;
LOW_PASS_CUT_OFF    = 10;

VELOCITY_THRESHOLD  = 5;
THRESHOLD_TIME      = FS_DS / 2;

% Step I: Low-pass filter
filteredSignal      = lowPassStream(signal,   FS, LOW_PASS_CUT_OFF);

% Step II: Downsample
downsampledSignal   = resample(filteredSignal, 1, 40);


% Step III: Velocity
velocity            = diff(downsampledSignal)./(1/FS_DS);

% Step IV: Threshold velocity
[onsets, offsets]   = consecEvents(abs(velocity - mean(velocity)) > VELOCITY_THRESHOLD, 1);
[onsets, offsets]   = pruneConsecEvents(onsets, offsets, THRESHOLD_TIME);

r = sum(offsets - onsets) / length(downsampledSignal) * 100;

% Plot Data
if nargout < 1
    figure;
    plot(1/FS_DS:1/FS_DS:length(downsampledSignal)/FS_DS, downsampledSignal);
    hold on;
    stem(onsets/FS_DS, 20.*ones(size(onsets)), 'ro');
    stem(offsets/FS_DS, 20.*ones(size(offsets)), 'ko');
    hold off;
    title(sprintf('%s; r: %0.3f', titleText, r));
    xlabel('Time (s)');
    ylabel('Velocity (degrees/sec)');
end
end
