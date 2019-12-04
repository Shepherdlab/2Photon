function [xyczt_img, metadata] = import_PrairieTif(varargin)
% SYNTAX: [xyczt_img metadata] = import_PrairieTif(varargin)
%   import_PrairieTif('C:\Images\SingleScan001-1')
%
% METHOD: Imports an image taken from a prairie CLSM into MATLAB memory
% with an axis convention of xyczt.  By default, removes any singleton 
% dimensions (see ALGORITHM). Also parses some metadata from the xml file.
%
% INPUT: 
% img_full_path [o, string]: path to the folder containing prairie tif 
%   files, xml file, and config file.  If not provided, prompts user to 
%   browse.
% channel_index [p, 1xnchan numeric]: index which reorder channels, the
%   value of each element is the index for each channel.  Extra channels 
%   are zero padded (i.e. can be used to force 2 channel image to be RGB by
%   specifying [1 2 3]).  Defaults to not reordering index or padding to
%   RGB.
% DataType [p, class]: the image will be scaled and casted to this
%   class denoted by this string. Default: uint16.
% SQUEEZE_DIMS [p, bool]: removes singleton dimensions of imported image.  
%   Default = true.
%
% OUTPUT:
% xyczt_img [img]: imported image, singleton dimensions may be
%   removed.
% metadata [struct]: a MATLAB struct containing a subset of the metadata
%   found in the xml file.  These are gleaned using regular expressions.  
%   See below for explanation.
%
% EXAMPLE:
%  [xyczt_img, metadata] = import_PrairieTif('C:\Images\ZSerie-015');
%  [xyczt_img, metadata] = import_PrairieTif();
%
% ALGORITHM: 
%   The img metadata is parsed from the xml file included with each prairie
% image.  Instead of parsing the entire xml file (which takes WAY TOO LONG 
% in MATLAB, even for small xml files), relevant metadata is gleaned using 
% regular expressions in token mode.  Code can easily be added to get more 
% metadata if desired.  
%   After loading the metadata, the individual prairie tif files are loaded
% into a 5D image of the form xyczt.  For our system, Prairie images are 12
% bit.  Since there is no native 12 bit dataype in matlab the images are 
% cast and scaled to 16 bit.  This behavior can be changed with the INPUT
% parameters.
%   By default (SQUEEZE_DIMS), singleton dimensions of the image are 
% removed, so a single image would be xyc, a z stack would be xycz, a time
% series z stack xyczt, and a time series single image would be xyct. 
% Disabling this will return a xyczt image regardless of whether those 
% dimensions are singleton or not.
%
% ASSUMPTIONS:
% The code assumes the bit depth and pixel dimensions remain constant
% across images.
%
% NOTES:
%  MATLAB Regular Expressions: mostly follows the standard regular expression
% syntax.  One exception found so far is '.' matches to any character,
% including newline.  So '.*' goes ACROSS newline characters, which is not
% default behavior. To make '.*' behave in correct way the easiest
% modification is '.*?' which means 'repeat matching of any char in a 
% RELUCTANT fashion ACROSS newline characters'.  So it matches to a single
% line if possible.
% 
%  Fields of 'metadata' structure.
%   .x_stg: x axis stage coordinate in um
%   .y_stg: y axis stage coordinate in um
%   .middle_z_stg: z axis stage coordinate of the middle z slice
%   .dx_um: x axis um/pixel ratio
%   .dy_um: y axis um/pixel ratio
%   .dz_um: um distance between z slices in um, NAN if single image.
%
%
% DEPENDENCIES: none
%   Tested MATLAB R2011a win32
%
% AUTHOR: Bruce A. Corliss, Scientist III
%
% LICENSE: BSD-3.  See Below.
%
% ORGANIZATION: GrassRoots Biotechnology
%
% CONTACT:
%  bruce.a.corliss@gmail.com
%  bruce.corliss@grassrootsbiotech.com
%
% DATE: 6/15/2011
%
% VERSION: 1.0.1

% LICENSE
% Copyright (c) 2011, GrassRoots Biotechnology
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

%% ARGUMENT PARSING
p = inputParser;
p.addOptional('img_full_path', '',@ischar);
p.addParamValue('channel_index', [], @isnumeric);
p.addParamValue('DataType', 'uint16', @char);
p.addParamValue('SQUEEZE_DIMS', 1, @(x)(x==0 || x==1));
p.StructExpand = false;  
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

% Prompt for manual browse to image folder if null input
if isempty(img_full_path); 
    img_full_path  = uigetdir('Select Folder containing Prairie Image');
    if(img_full_path ==0); return; end
end


%% If image folder or tif or XML files DNE, return empty
if isempty(dir(img_full_path)) || numel(dir([img_full_path '/*.tif']))==0 ...
        || numel(dir([img_full_path '/*.xml']))==0; 
    xyczt_img = []; return;
end


%% Make path *nix compatible, extract name of image
img_full_path = regexprep(img_full_path, '\', '/');
[~, img_name] = fileparts(img_full_path);


%% Read in metadata from xml file
text = fileread([img_full_path '/' img_name '.xml']);


%% Tile stitching is not supported yet.  Error if detected
TILE = strcmp(regexp(text, 'xYStageGridDefined="(\w*)"', 'tokens','once'), 'True');
if TILE; errordlg(['Tile stitching not supported yet, single composite '...
        'image will be constructed.']); end


%% Find dimensions of 5D image
% XY dimensions: Entries assumed to be the same for pixel dimensions
xdim = str2double(regexp(text, ...
    'Key key="pixelsPerLine".*?value="(\d*)"', 'tokens','once'));
ydim = str2double(regexp(text, ...
    'Key key="linesPerFrame".*?value="(\d*)"', 'tokens','once'));

% Z axis dimension
z_cell =  regexp(text, 'Frame relative.*?index="(\d*)"', 'tokens');
zslice = cellfun(@(x) str2double(x{1}),z_cell);

% T axis dimension
t_cell =  regexp(text, 'Sequence type=.*?cycle="(\d*)"', 'tokens');
tpoints = cellfun(@(x) str2double(x{1}),t_cell);
if tpoints == 0; tpoints=1; end

% Channel dimension (img names found on same line, parsed also)
ch_img_names_cell =  regexp(text, ...
    'File channel="(\d*)".*?filename="(.*?)"', 'tokens');
ch = cellfun(@(x) str2double(x{1}),ch_img_names_cell);
img_names = cellfun(@(x) x{2},ch_img_names_cell, 'UniformOutput', 0);

% Bit depth of img in filesystem
tif_bit_depth = str2double(regexp(text, ...
    '<Key key="bitDepth".*?value="(\d*)"', 'tokens','once'));


%% Parse extra subset of metadata
% XYZ Coordinates of each img
% It is assumed that the xy stage coordinates are constant throughout the
% image (tiling is not supported yet).  Sometimes the data is not written
% correctly to the xml file and the value is listed as zero.  To deal with
% this the mode of the entries are taken.
metadata.x_stg = mode(cellfun(@(x) str2double(x{1}), regexp(text, ...
    '<Key key="positionCurrent_XAxis".*?value="([-\d\.]*?)"', 'tokens')));
metadata.y_stg = mode(cellfun(@(x) str2double(x{1}), regexp(text, ...
    '<Key key="positionCurrent_YAxis".*?value="([-\d\.]*?)"', 'tokens')));
metadata.middle_z_stg = median(unique(cellfun(@(x) str2double(x{1}), ...
    regexp(text, '<Key key="positionCurrent_ZAxis".*?value="([-\d\.]*?)"', 'tokens'))));
metadata.dx_um = mode(cellfun(@(x) str2double(x{1}), regexp(text, ...
    '<Key key="micronsPerPixel_XAxis".*?value="([-\d\.]*?)"', 'tokens'))); 
metadata.dy_um = mode(cellfun(@(x) str2double(x{1}), regexp(text, ...
    '<Key key="micronsPerPixel_YAxis".*?value="([-\d\.]*?)"', 'tokens'))); 
metadata.dz_um = mode(diff(cellfun(@(x) str2double(x{1}), regexp(text,...
    '<Key key="positionCurrent_ZAxis".*?value="([-\d\.]*?)"', 'tokens'))));


%% Create 1x1 mapping between img_name, channel, z slice, timepoint
% Channel index
ch_ind = ch;  
if ~isempty(channel_index)
    if ~isempty(setdiff(ch_ind, channel_index)); error(['Specified channel index must'...
            ' contain all channel index indices parsed from xml file']);
    end
    % Reassign channel indices
    switch_cell = arrayfun(@(x) x == ch_ind, 1:numel(unique(ch)), ...
        'UniformOutput', 0);
    for n = 1:numel(switch_cell); ch_ind(n) = channel_index(n); end
end

% Zslice and timpoint index need to be repeated to match img_names elements
flat = @(x) x(:);
z_ind = flat(repmat(zslice, [numel(unique(ch_ind)) 1]));
t_ind = flat(repmat(tpoints, [numel(unique(z_ind))*numel(unique(ch_ind)) 1]));


%% Clear xml text form memory, not needed for final img import
clear text;


%% Initialize image specified datatype
xyczt_img = zeros(xdim, ydim, max(horzcat(ch, channel_index)), ...
    max(zslice), max(tpoints), DataType);


%% Read individual tif files into 5D img
for n = 1:numel(img_names)
    xyczt_img(:,:,ch_ind(n), z_ind(n), t_ind(n)) = imread([img_full_path...
        '/' img_names{n}]) * double(intmax(DataType)/2^tif_bit_depth);
end


%% Remove extra singleton dimensions
if SQUEEZE_DIMS
    xyczt_img = squeeze(xyczt_img);
end


end