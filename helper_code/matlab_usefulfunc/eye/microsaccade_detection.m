function ms = microsaccade_detection(eye_x, eye_y, samplingRate, algorithm, fig)
%% 
% detect mscade events in the given x,y eye-positions
% during fixation by using a specified algorithm
%
% INPUT: eye_x ... vector of horizontal eye position in one trial of fixation
%        eye_y ... vector of vertical eye position in one trial of fixation
%        samplingRate ... sampling rate in the measurement system
%        algorithm ... 'Nienborg'; Nienborg & Cumming (2006) ' http://www.jneurosci.org/content/26/37/9567.long' 
%                  ... 'Engert'; Engbert & Kliegl (2003) 'http://www.sciencedirect.com/science/article/pii/S0042698903000841'
%                  ... 'Joshi'; Joshi et al (2016) 'http://www.sciencedirect.com/science/article/pii/S089662731501034X'
%                  ... 'hafed'; Hafed et al (2011) 'http://www.jneurosci.org/content/31/43/15219.long'
%        fig ... 1; validation plot
% OUTPUT: ms ... matlab structure file storing relevant infos
%
% NOTE: As a detection algorithm 'Engert' is recommended, as this is the
% one widely used.
%
% EXAMPLE: ms = mscade(eye_x, eye_y, 500, 'Nienborg', 1);
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++

%%
% deal with inputs
if nargin < 2 
    error('At least eye-positions for x and y must be given as the first and second input arguments.')
end
if nargin < 3; samplingRate = 500; end
if nargin < 4; algorithm = 'Nienborg'; end
if nargin < 5; fig = 1; end

%%
% match the length of eye_x and eye_y, if different
if length(eye_x) ~= length(eye_y)
  disp('eye-position x and y have different length.')
  % force eye_x and eye_y to be the same length
    while abs(length(eye_x)-length(eye_y))>0
        if length(eye_x) > length(eye_y)
            eye_x(end) = [];
        else
            eye_y(end) = [];
        end
    end
    disp('Length of eye-position x and y were matched.')
end
len_eye = length(eye_x);

%%
% define a detection algorithm
ms.algorithm.maxamp = 1; % deg
switch lower(algorithm)
    case  'nienborg' % Nienborg & Cumming (2006)
        ms.algorithm.name = 'Nienborg & Cumming (2006)';
        ms.algorithm.velocity = 12; % deg/s
        ms.algorithm.smooth = 0;
        ms.algorithm.acceleration = nan;
        ms.algorithm.minduration = nan;

      case 'engbert' % Engbert & Kliegl (2003); Engbert & Mergenthaler (2006)
        % transform eye-positions into velocities
        vel_x = zeros(1,len_eye);
        vel_y = zeros(1,len_eye);
        for i = 3:len_eye-2
              vel_x(i) = samplingRate*(eye_x(i+2)+eye_x(i+1)-eye_x(i-1)-eye_x(i-2))/6;
              vel_y(i) = samplingRate*(eye_y(i+2)+eye_y(i+1)-eye_y(i-1)-eye_y(i-2))/6;
        end

        % SD of the velocity distribution as the detection threshold
        sigma = zeros(2,1);
        sigma(1) = median(vel_x.^2) - (median(vel_x))^2;
        sigma(2) = median(vel_y.^2) - (median(vel_y))^2;
        gamma = 6; % default used in the original paper was 5 or 6
        thre_x = gamma*sigma(1);
        thre_y = gamma*sigma(2);
        
        % store algorithm info
        ms.algorithm.name = 'Engbert & Kliegl (2003)';
        ms.algorithm.velocity = [thre_x, thre_y]; % deg/s
        ms.algorithm.smooth = 1;
        ms.algorithm.acceleration = nan;
        ms.algorithm.minduration = 6; % ms
            
    case 'joshi' % Joshi et al., 2016
        ms.algorithm.name = 'Joshi et al (2016)';
        ms.algorithm.velocity = 15; % deg/s
        ms.algorithm.smooth = 0;
        ms.algorithm.acceleration = nan;
        ms.algorithm.minduration = 6; % ms
        
   case 'hafed' % Hafed et al., 2011        
        ms.algorithm.name = 'Hafed et al (2011)';
        ms.algorithm.velocity = 8; % deg/s
        ms.algorithm.smooth = 0;
        ms.algorithm.acceleration = 550;
        ms.algorithm.minduration = nan; % ms
end

%%
% detection ms by velocity
event = zeros(1,len_eye);
if ms.algorithm.smooth == 0
    % compute instantaneous eye velocity
    vel_x = [0 diff(eye_x)*samplingRate];
    vel_y = [0 diff(eye_y)*samplingRate];
    vel = sqrt(vel_x.^2 + vel_y.^2);
    event(vel > ms.algorithm.velocity) = 1;      
elseif ms.algorithm.smooth == 1
    % detect ms events based on the ellipse equation          
    event((vel_x.^2)/(ms.algorithm.velocity(1)^2) + ...
        (vel_y.^2)/(ms.algorithm.velocity(2)^2) > 1) = 1;
end

% detection ms by acceleration
if ~isnan(ms.algorithm.acceleration)
    vel = sqrt(vel_x.^2 + vel_y.^2);
    acc = [0 diff(vel)];
    event(acc > ms.algorithm.acceleration) = 1;
end

% detection ms by minimum duration
minsamples = ceil(ms.algorithm.minduration/((1/samplingRate)*1000));
if isnan(minsamples)
    minsamples = 2;
end

% apply minimum duration criteria
[sacc_start, sacc_end, count] = consecutive_ones(event);
sacc_start(count < minsamples) = [];
sacc_end(count < minsamples) = [];
count(count < minsamples) = [];

% store eye-data into a structure
if isempty(count)
    ms.counts = 0;
else
    ms.counts = length(sacc_start);
end
if ms.counts > 0
    ms.amp = zeros(1, ms.counts);
    ms.peakv = zeros(1, ms.counts);
    ms.duration = zeros(1, ms.counts);
    ms.angle = zeros(1, ms.counts);
    for i = 1:ms.counts
        amp = [];
        ang = [];
        j = 0;
        while sacc_start(i) + j < sacc_end(i)
            amp = [amp, sqrt((eye_x(sacc_start(i)+j+1) - eye_x(sacc_start(i)+j)).^2 ...
            + (eye_y(sacc_start(i)+j+1) - eye_y(sacc_start(i)+j)).^2)];
            ang = [ang, atan2(eye_y(sacc_start(i)+j+1) - eye_y(sacc_start(i)+j), ...
            eye_x(sacc_start(i)+j+1) - eye_x(sacc_start(i)+j))*180/pi];
            j = j + 1;
        end
        ms.amp(i) = sum(amp);        
        ms.angle(i) = median(ang);
%         ms.amp(i) = sqrt((eye_x(sacc_end(i)) - eye_x(sacc_start(i))).^2 ...
%             + (eye_y(sacc_end(i)) - eye_y(sacc_start(i))).^2);
%         ms.angle(i) = atan2(eye_y(sacc_end(i)) - eye_y(sacc_start(i)), ...
%             eye_x(sacc_end(i)) - eye_x(sacc_start(i)))*180/pi;
        ms.peakv(i) = max([sqrt(vel_x(sacc_start(i):sacc_end(i)).^2 ...
            + vel_y(sacc_start(i):sacc_end(i)).^2)]);
        ms.duration(i) = (count(i))/samplingRate;
    end
    
%     % maximum amplitude criteria
%     idx = ms.amp <= ms.algorithm.maxamp;
%     ms.counts = sum(idx==1);
%     fieldnames = {'amp', 'peakv', 'duration', 'angle'};
%     for i = 1:length(fieldnames)
%         ms.(fieldnames{i}) = ms.(fieldnames{i})(idx);
%     end
else
    ms.amp = nan;
    ms.peakv = nan;
    ms.duration = nan;
    ms.angle = nan;
end
ms.velocity_x = vel_x;
ms.velocity_y = vel_y;
ms.event = event;

% validate with plot
if fig==1
    if ms.counts == 0
        disp('no saccade is detected.')
        return
    end
      figure;
      lw = 0.25;
      % time domaincsa
      subplot(4,6,1:6)
      time = [1:length(eye_x)]/samplingRate;
      plot(time, eye_x, '-','Color',[0.5 0.5 0.5])
      ylabel('x (deg)')
      set(gca,'box','off'); set(gca,'TickDir','out');  
      subplot(4,6,7:12)
      plot(time, eye_y, '-','Color',[0.5 0.5 0.5])
      ylabel('y (deg)')
      xlabel('time (sec)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      % position space
      subplot(4,6,13:15)
      plot(eye_x,eye_y,'-','Color',[0.5 0.5 0.5], 'linewidth', lw)
      xlabel('horizontal position (deg)')
      ylabel('vertical position (deg)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      % velocity space
      subplot(4,6,16:18)
      plot(vel_x, vel_y, '-','Color',[0.5 0.5 0.5], 'linewidth', lw)
      hold on;
      if ms.algorithm.smooth == 1
          draw_ellipse(thre_x, thre_y, [0 0])
      else
          th = 0:pi/50:2*pi;
          xunit = ms.algorithm.velocity * cos(th);
          yunit = ms.algorithm.velocity * sin(th);
          plot(xunit, yunit,'--k');
      end
      xlabel('horizontal velocity (deg/s)')
      ylabel('vertical velocity (deg/s)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      % draw saccadic traces
      map = lines(ms.counts);
      for i = 1:ms.counts
          k = 0;
          while sacc_start(i) + k < sacc_end(i)
              subplot(4,6,1:6)
              hold on;
              plot(time(sacc_start(i)+k:sacc_start(i)+k+1), eye_x(sacc_start(i)+k:sacc_start(i)+k+1),...
                  '-','Color',map(i,:))
              subplot(4,6,7:12)
              hold on;
              plot(time(sacc_start(i)+k:sacc_start(i)+k+1), eye_y(sacc_start(i)+k:sacc_start(i)+k+1),...
                  '-','Color',map(i,:))
              subplot(4,6,13:15)
              hold on;
              plot([eye_x(sacc_start(i)+k) eye_x(sacc_start(i)+k+1)], [eye_y(sacc_start(i)+k) eye_y(sacc_start(i)+k+1)],...
                      '-','Color',map(i,:), 'linewidth', lw)
              subplot(4,6,16:18)
              hold on;
              plot([vel_x(sacc_start(i)+k) vel_x(sacc_start(i)+k+1)], [vel_y(sacc_start(i)+k) vel_y(sacc_start(i)+k+1)],...
                      '-','Color',map(i,:), 'linewidth', lw)
              k = k + 1;
          end
      end

      subplot(4,6,[19 20])
      plot(ms.amp, ms.peakv,'.k')
      set(gca, 'XScale', 'log')
      set(gca, 'YScale', 'log')
      xlabel('amplitude (deg)')
      ylabel('peak velocity (deg/s)')
      set(gca,'box','off'); set(gca,'TickDir','out');  

      subplot(4,6,[21 22])
      polarhistogram(ms.angle*pi/180)
      set(gca,'box','off'); set(gca,'TickDir','out'); 
      
      subplot(4,6,[23 24])
      histogram(1000*ms.duration')
      xlabel('duration (ms)')
      set(gca,'box','off'); set(gca,'TickDir','out'); 
end

% subfunction
function [start_idx, end_idx, count] = consecutive_ones(vector)
% https://de.mathworks.com/matlabcentral/fileexchange/34914-consecutive-ones
%[Function Description]
%Finds the number of consecutive ones in a binary signal. Returns the
%starting and ending points of the consecutive set and also the count. Since it
%does not use any for loop. It is pretty fast
%
%[Input]
%vector = A binary signal
%[Output]
%
%Author - Shreyes

temp = diff([0 vector 0]);
start_idx = find(temp == 1);
end_idx = find(temp == -1);
count = end_idx - start_idx;
end_idx = end_idx -1;

function draw_ellipse(h, v, cc)
t=-pi:0.01:pi;
x = cc(1) + h*cos(t);
y = cc(2) + v*sin(t);
plot(x,y, '--k')