
% CHECK_SPEC  Checks if a configuration files exists for given spec and loads.
%  
%  AUTHOR: Timothy Sipkens, 2024-01-08

function prop = check_spec(spec)

% Get folder where this script is located. 
% Allows for running this script when added as a submodule.
fd = fileparts(mfilename('fullpath'));

% Get corresponding file name for configuration file.
fn = [fd, filesep, 'config', filesep, spec, '.json'];  % file name of potential config

if isfile(fn)
    prop = load_config([fd, filesep, 'config', filesep, spec, '.json']);
else
    prop = -1;  % return -1 if file not found
end

end


%=========================================================================%
% LOAD_CONFIG  Loads or update settings from configuration file (JSON). 
%  Files are loaded in order supplied, overwriting properties where
%  relevant (e.g., in selecting the setting for k-means segmentation).
function prop = load_config(fnames, prop)

if isempty(fnames); fnames = {}; end
if ~iscell(fnames); fnames = {fnames}; end

if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = struct(); end

% If default prop is given as a file name.
% Load this file first. 
if isa(prop, 'char')
    prop = load_files({prop}, struct());
end

% If fnames is not empty, load the second file. 
if ~isempty(fnames)
    prop = load_files(fnames, prop);
end

end



%== LOAD_FILES ==========================================%
%   Loading settings from configuration file (JSON). 
function prop = load_files(fnames, prop)

for ii=1:length(fnames)
    
    prop0 = read_json(fnames{ii});  % read new settings
    
    f = fieldnames(prop0);
    
    % Copy (or overwrite) existing settings.
    for jj = 1:length(f)
        
        % Attempt to interpret Matlab expressions.
        prop0.(f{jj}) = interpret(prop0.(f{jj}));
        
        % Overwrite existing configuration.
        prop.(f{jj}) = prop0.(f{jj});
    end
    
end

end


%== INTERPRET ===========================================%
%   Attempt to interpret Matlab expressions.
function e0 = interpret(f)

if isa(f, 'char')
    [e0, success] = str2num(f);
    if success; return; end

    try  % try 'eval(...)'
        e0 = eval(f);
    catch
        e0 = f;
        return;  % just continue
    end
else
    e0 = f;
end

end


%== READ_JSON ==========================================%
%   Read JSON structured configuration files. 
%   Allows for C++ or Javscript style commenting.
%  
%   AUTHOR: Timothy Sipkens, 2021-04-20

function results = read_json(file)

fid = fopen(file);
raw = fread(fid, inf);  % raw file contents
str = char(raw');  % transform to char
fclose(fid);

% Remove comments.
str = erase(erase(eraseBetween( ...
    erase(eraseBetween(str, "//", newline), "//"), ...
    "/*", "*/"), "/*"), "*/");

results = jsondecode(str);

end



