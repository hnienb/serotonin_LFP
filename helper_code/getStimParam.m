function [param, vals] = getStimParam( ex )
% getStimParam returns the one stimulus that was modulated during the
% experiment, i.e. that is in the ex.Trials file
% 
% @CL

param = ex.exp.e1.type;    %<- stimulus dimension (contrast, orientation, ...)
vals = unique([ex.Trials.(param)]); %<- stimulus samples

vals = sort(vals); % sort the parameters in ascending order

if isfield(ex.exp.e1, 'blank')
    if any(vals == ex.exp.e1.blank)
        vals = [vals(end) vals(1:end-1)]; % make blank first    
    end
end

