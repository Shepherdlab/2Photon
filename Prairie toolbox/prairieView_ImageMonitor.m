function img_list = prairieView_ImageMonitor(xyz_loc, out_path, varargin)
% SYNTAX: img_list = prairieView_ImageMonitor(xyz_loc, out_path, varargin)
%
% METHOD: Given a list of locations (xyz_loc), acquires a zstack of at each
% with the configuration specified (config_name) that is saved as a file on
% the hard disk (prairie_config_path). Images are saved to hard disk
% (out_path), with the img_base_name and file_iternum specified.
% Custom z stack settings can be specified (zstack_nzdz), or the current
% settings are used.
% 
% NOTE: For this code to "just work" as intended, a couple of controls have
% to be set in PrairieView
%   1. In the Z Series Tab: Check [ ] At all XY stage locations
%       If this is not checked prairie view will treat the xyz locations to be
%       the start of the z stack instead of the middle, even if this is
%       specified in the settings.
%   2. Set Start Position to [ 0 ]
%   3. Uncheck [x] Adjust PMT & Laser
%   4. Select ( ) Define in the Stop Position Section so it changes to 
%       (*) Calculate
%   5. Other settings are default.
%
% INPUT:
%   xyz_loc [r, nx3 numeric]: xyz coordinates of each location to be imaged.
%   out_path [r, string]: path where images are to be saved.
%   config_name [r, string]: name of config file to be loaded (this is
%       assumed to be saved to disk (in prairie_config_path) as a file with
%       the same name.
%   file_iternum [p, int]: filename iteration number for first image.
%   prairie_config_path [p, str]: path to where imaging configs are stored
%       as files.
%   img_base_name [p, str]: base file name of images saved by prairie view.
%   zstack_nzdz [p, 1x2 numeric]: two element numeric vector, use to
%       specify the parameters for a z stack: number of z slices and um 
%       spacing between z slices. This info is not saved with imaging
%       configurations.
%   img2scope_axdir [p, 1x3 numeric]: this xyz vector is multiplied to the 
%       xyz coordinates specified. This is used to switch the sign of the
%       coordinates if the prairie view stage has different directions for
%       xyz compared to the one used for calculations.
%   CHECK_OUT_DIR [p, bool]: true makes program check output directory if
%       images with the same base name (and right # of images) exist. If
%       they do, then imaging is skipped and the filenames are returned. 
%   HIDE_WAITBAR [p, bool]: true hides a waitbar that appears to show
%       imaging progress.
% 
% OUTPUT:
%   img_list [cell]: list of image names.
%
% DEPENDENCIES: 
%   prairieView_Com
%   
% AUTHOR: Bruce Corliss, Scientist III
% DATE: 6/15/2012

% LICENSE
% Copyright (c) 2012, GrassRoots Biotechnology
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
%     Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     Neither the name of the GrassRoots Biotechnology nor the names of its 
%       contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
% AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
% USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if getappdata(0, 'VERBOSE'); fprintf('%s\n', mfilename); end

%% PARAMTERS
% Enable pause function
pause on
digits  = @(x) length(sprintf('%.f', abs(x)));

%% ARGUMENT PARSING
p = inputParser;
p.addRequired('xyz_loc', @(x) size(x,2)==3 && isnumeric(x));
p.addRequired('out_path', @ischar);
p.addParamValue('config_name', 0, @(x) ischar(x));
p.addParamValue('file_iternum', 1, @(x) isnumeric(x));
p.addParamValue('prairie_config_path', '', @(x) ischar(x));
p.addParamValue('img_base_name', 'pv_img', @(x) ischar(x));
p.addParamValue('zstack_nzdz', [], @(x) (isnumeric(x) && numel(x)==2) || isempty(x));
p.addParamValue('img2scope_axdir', [1 1 1], @(x) (isnumeric(x) && numel(x)==3));
p.addParamValue('CHECK_OUT_DIR', 1, @(x)(x==0 || x==1));
p.addParamValue('HIDE_WAITBAR', 0, @(x)(x==0 || x==1));
p.addParamValue('DEBUG_FLAG', 0, @(x)(x==0 || x==1));
p.StructExpand = false;  
p.parse(xyz_loc, out_path, varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

% out_path required input
if isempty(out_path); error([mfilename ': CheckDir required']); end


%% Print event type and number
nim = size(xyz_loc,1);
fprintf('\t- PV Monitor Event: %d images.\n',nim); 

% Check to see if img folders have been moved to out_path
if CHECK_OUT_DIR
    if numel(dir([out_path '/' img_base_name '*'])) == size(xyz_loc,1);
        
        fprintf('\t- Images acquired already.\n');
        img_list = arrayfun(@(x) x.name, dir([out_path '/' img_base_name '*']), 'UniformOutput', 0);
        return;
    end
end


%% Load config File, fit image to window
config_path = [prairie_config_path '/' config_name '.cfg'];
prairieView_Com(['-lcf "' config_path '" -iwf']);


%% Set base name for saved images, file iternation # in base name, save path
prairieView_Com(['-fn Zseries "' img_base_name '"'...
    ' -fi Zseries ' sprintf('%u', file_iternum)...
    ' -p "' out_path '"']);


%% Set z stack paramters
if ~isempty(zstack_nzdz)
    % Z stack slices
    prairieView_Com(['-zsn ' sprintf('%d', (zstack_nzdz(1)))...
        ' -zsz ' sprintf('%d', (img2scope_axdir(3)*zstack_nzdz(2)))]);
end


%% Go through all location in list (or manually specify coordinates)
tb = repmat(' ', 1, 5);
if ~HIDE_WAITBAR
        hw = waitbar(0, sprintf(['%s\n%s' tb '%s' tb '%s  %s'], regexprep(config_name, '_',' '),...
            'ETA:',' -- : -- : --', num2str(nim), 'Images'), 'Name', ...
            'Microscope Acquiring Images','WindowStyle', 'modal');
end

% Initiate vector to record system time after each image for estimation of
% completion
datenum_vect = zeros(1,nim+1);
datenum_vect(1) = datenum(clock);

fprintf('\t- Imaging %u locations with %s.\n', nim, config_name);
for n = 1:nim

    % Clear stage positions, add new one, move to it
    prairieView_Com([...
        '-fi Zseries ' sprintf('%u', file_iternum+n-1) ' -spc ' ...
        ' -spa ' num2str(xyz_loc(n,:)) ...
        ' -mtsp ' sprintf('%u',1)]);
    
    % Take Z stack, wait for completion
    prairieView_Com('-zs -w');
    
    %% Update waitbar
    datenum_vect(n+1)= datenum(clock); 
    % Fit a linear equation for eta
    pfit = polyfit(1:n+1, datenum_vect(1:n+1),1);
    
    % Get eta
    eta = round(datevec((pfit(1)*nim+1 + pfit(2)) - max(datenum_vect)));
    
    % Pretty print numbers to string for display
    eta_cell = arrayfun(@(x) sprintf('%2.2u', x), eta, 'UniformOutput', 0);
    eta_str = [eta_cell{4} ' : ' eta_cell{5} ' : ' eta_cell{6}];
    nim_str = ...
        sprintf(['%' num2str(digits(nim)) '.u / %' num2str(digits(nim)) '.u'],...
        n, nim);
    
    if ~HIDE_WAITBAR
            % Update waitbar to show eta
            waitbar(n/nim, hw, sprintf(['%s\n%s' tb '%s' tb tb tb '%s %s'],...
                regexprep(config_name, '_',' '), 'ETA:', eta_str, 'Imaged:',nim_str));
    end
end;

if ~HIDE_WAITBAR; close(hw); end


%% Clear stage positions in prairie view, Set no action after scan
prairieView_Com(['-spc' ' -as']);


% Obtain list of files (could be written smarter)
img_iter_range = file_iternum:file_iternum+nim-1;
item_list = arrayfun(@(x) x.name, dir([out_path '/' img_base_name '*']), 'UniformOutput', 0);
isdir_bv = arrayfun(@(x) x.isdir, dir([out_path '/' img_base_name '*']));
img_dir_list = item_list(isdir_bv);
img_iter_observed = cellfun(@(x) str2double(regexp(x, '.*-(\d)*', 'tokens', 'once')), img_dir_list);
img_list = img_dir_list(intersect(img_iter_range,img_iter_observed));


% keyboard
end



