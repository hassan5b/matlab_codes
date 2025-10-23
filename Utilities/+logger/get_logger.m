%% A function used to get access to a single common instance of a logger.

% Authors:  
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 21st of July 2021

function [obj, delete_logger_handle] = get_logger()
  
    persistent l;
    persistent owner_count;
    if isempty(l)
        l = logger.logger();
    end
    obj = l;
    owner_count = owner_count + 1;
    
    delete_logger_handle = @() delete_logger();
    
    function delete_logger()
        owner_count = owner_count - 1;
        if owner_count == 0
            delete(obj);
            l = [];
            clear('obj');
        end
    end
end

