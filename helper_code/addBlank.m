% 
% add...
% - ...blank value 
% - ...tolerance information, if entries do not exactly match the set value
%   settings
% 
% @author Corinna Lorenz
% @created 22.09.2015
% @lastchanged 22.09.2015


% all occuring stimulus parameters
allvals = unique([ex.Trials.(ex.exp.e1.type)]); 

% add the value for blanks
ex.exp.e1.blank = unique([ex.Trials([ex.Trials.st]==0).(ex.exp.e1.type)]);
ex.exp.e1.blankname = 'blank';


% add blank and tolerance
if length(ex.exp.e1.blank) > 1
    ex.exp.e1.blank = ex.exp.e1.blank(1);
    [ex.Trials([ex.Trials.st]==0).(ex.exp.e1.type)] = deal(ex.exp.e1.blank);
elseif isempty(ex.exp.e1.blank) 
    % add e1 blank in case there was none
     ex.exp.e1.blank = 100001;     
end

% add e2 dummy
if isfield(ex, 'exp') && ~isfield(ex.exp, 'e2')
    ex.exp.e2.type  = 'me';
    ex.exp.e2.vals  = unique([ex.Trials.me]);
end

clearvars allvals;