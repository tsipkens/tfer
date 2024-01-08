
% LOAD_SPEC  Checks if a configuration files exists for given spec and loads.
%  
%  AUTHOR: Timothy Sipkens, 2024-01-08

function prop = load_spec(spec)

% Get folder where this script is located. 
% Allows for running this script when added as a submodule.
fd = fileparts(mfilename('fullpath'));

% Get corresponding file name for configuration file.
fn = [fd, filesep, 'config', filesep, spec, '.yaml'];  % file name of potential config
fnl = [fd, filesep, 'config_local', filesep, spec, '.yaml'];  % allow for local config

if isfile(fnl)  %try local load first
    prop = read_yaml(fnl);
elseif isfile(fn)
    prop = read_yaml(fn);

else
    prop = -1;  % return -1 if file not found
end

end



% READ_YAML  Simple utility to read simple YAML files.
%  
%  AUTHOR: Lloyd Russell, 2017
%  REPO: https://github.com/llerussell/ReadYAML
%  MODIFIED: Timothy Sipkens, 2021-04-19
function prop = read_yaml(file_path, prop)

% Parse inputs (currently never invoked).
if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = struct(); end

% If file of presets does not exist.
if and(~isfile(file_path), ~isfile(['prop', filesep, file_path]))
    error(['PMA properties not available for provided case. ',...
            'Try a different string in the `prop_pma(str)` call.']);
end


% Read file line by line.
fid = fopen(file_path, 'r');
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);

% Remove empty lines.
data = deblank(data{1});
data(cellfun('isempty', data)) = [];

% Prepare final results structure.
results = [];

% Parse the contents (line by line).
for i = 1:numel(data)
    
    % extract this line
    thisLine = data{i};
    
    % ignore if this line is a comment
    if strcmpi(thisLine(1), '#')
        continue
    end
    
    % find the seperator between key and value
    sepIndex = find(thisLine==':', 1, 'first');
    
    % get the key name (remove whitespace)
    key = strtrim(thisLine(1:sepIndex-1));  
    
    % get the value, ignoring any comments (remove whitespace)
    value = strsplit(thisLine(sepIndex+1:end), '#');
    value = strtrim(value{1}); 
    
    % attempt to convert value to numeric type
    [convertedValue, success] = str2num(value);
    if success
        value = convertedValue;
    end
    
    % store the key and value in the results
    results.(key) = value;
end

% Append to supplied structure.
% Overwrite duplicate fields. 
f_results = fields(results);
for ii = 1:length(f_results)
    prop.(f_results{ii}) = results.(f_results{ii});
end

end


