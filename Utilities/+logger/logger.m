%% A simple logger which supports various levels of logging and prints to the command prompt
% Highly trimmed and modified version of 
% https://github.com/optimizers/logging4matlab

% Authors:  
%   Neelay Junnarkar, <neelay.junnarkar -AT- berkeley.edu>, EECS, UC Berkeley
%   Pierre-Jean Meyer, <pierre-jean.meyer -AT- univ-eiffel.fr>, COSYS-ESTAS, Univ Gustave Eiffel
% Date: 2nd of December 2021

classdef logger < handle
   
    properties (Constant)
        % Increasing number implies increasing priority.
        % For a given log level, all messages with priority at least
        % as large will be printed.
        ALL     = 0;
        RUNTIME = 1;
        INFO    = 2;
        WARNING = 3;
        ERROR   = 4;
        OFF     = 5;
        
        level_to_string = containers.Map(...
            {logger.logger.ALL,  logger.logger.RUNTIME,...
            logger.logger.INFO,  logger.logger.WARNING,...
            logger.logger.ERROR, logger.logger.OFF},...
            {'ALL', 'RUNTIME', 'INFO', 'WARNING', 'ERROR', 'OFF'}); 
        
        
        % File ID standard output (command window).
        STANDARD_OUTPUT = 1;
        
        % Date string format to prepend messages.
        % Year-month-day hour:minute:second
        date_string_format = 'yyyy-mm-dd HH:MM:SS';
    end
    
    properties (SetAccess = protected)
        log_level = logger.logger.ALL;
        file_id = logger.logger.STANDARD_OUTPUT;
    end
    
    methods
        
        % log_level should be one of the constants ALL, RUNTIME, INFO,
        % WARNING, ERROR, or OFF.
        % For a given log level, all messages of equal or higher log levels
        % will be displayed.
        function set_log_level(self, log_level)
            assert(ismember(log_level, [self.ALL, self.RUNTIME, self.INFO,...
                self.WARNING, self.ERROR, self.OFF]));
            self.log_level = log_level;
        end
        
        % Set log_location to:
        %   NaN: print to command window.
        %   file name: print to file (by appending if it already exists).
        %       e.g. log.set_log_location("logfile.log")
        function set_log_location(self, log_location)
            if ~isstring(log_location)
                self.file_id = self.STANDARD_OUTPUT;
                return;
            end
            
            self.file_id = fopen(log_location, 'a');
            
            % Default to standard output if fail to open file.
            if self.file_id == -1
                self.file_id = self.STANDARD_OUTPUT;
            end
        end
            
        
        % Log with log level RUNTIME.
        function runtime(self, format, varargin)
            self.log(self.RUNTIME, format, varargin{:});
        end
        
        % Log with log level INFO.
        function info(self, format, varargin)
            self.log(self.INFO, format, varargin{:});
        end
        
        % Log with log level WARNING.
        function warning(self, format, varargin)
            self.log(self.WARNING, format, varargin{:});
        end
        
        % Log with log level ERROR.
        function error(self, format, varargin)
            self.log(self.ERROR, format, varargin{:});
        end
    end
    
    methods (Access = protected)
        function log(self, log_level, format, varargin)
            if log_level < self.log_level
                return;
            end
            
            line = sprintf(format, varargin{:});
            
            if ~(log_level == self.ERROR && self.file_id == self.STANDARD_OUTPUT)              
                fprintf(self.file_id, '[%s] %s:\t%s',...
                    datestr(now(), self.date_string_format), self.level_to_string(log_level), line);
            end
            if log_level == self.ERROR
            	error(format, varargin{:});
            end
        end
    end
end

