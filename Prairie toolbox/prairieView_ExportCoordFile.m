function EXPORT_SUCCESS = prairieView_ExportCoordFile(coord_path, xyz_loc)
% SYNTAX: EXPORT_SUCCESS = prairieView_ExportCoordFile(coord_path, xyz_loc)
%
% METHOD: exports xyz coordinates to filesystem in format that can be
% loaded by prairie view.
%
% INPUT:
%   coord_path [r, path]: string of full file path where the coordinates
%       are to be saved.
%   xyz_loc [r, nx3 double]: for (n) locations, nx3 xyz list of stage
%       locations.
% 
% OUTPUT:
%   EXPORT_SUCCESS [boolean]: true if file successfully written to
%       filesystem
%
% DEPENDENCIES: nones
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

% Assume success
EXPORT_SUCCESS = 1;

% Print to stdout
fprintf('\t%s:\t%s\n', mfilename, coord_path);

% Create output directory id DNE
if isempty(dir(fileparts(coord_path))); mkdir(fileparts(coord_path)); end

% Attempt to open file for writing
fid = fopen(coord_path, 'w');
if fid == -1;EXPORT_SUCCESS = 0; fprintf('\t- Failed to write coord file.\n');
    return; end

% Write XML
fprintf(fid, '%s\n','<?xml version="1.0" encoding="utf-8"?>');
fprintf(fid, '%s\n','<StageLocations>');

CoordStr =@(n,x,y,z) ['<StageLocation index="' sprintf('%u', n) ...
    '" x="' sprintf('%.2f', x) '" y="' sprintf('%.2f', y) '" z="' ...
    sprintf('%.2f', z) '" />'];

% Write each location in xml format
for n=1:size(xyz_loc,1)
   fprintf(fid, '%s\n', CoordStr(n, xyz_loc(n,1),xyz_loc(n,2),xyz_loc(n,3))); 
end
fprintf(fid, '</StageLocations>');

% Close file
fclose(fid);

end