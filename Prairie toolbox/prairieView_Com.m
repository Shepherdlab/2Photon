function status = prairieView_Com(cmd_str, varargin)
% SYNTAX: status = prairieView_Com(cmd_str, varargin)
%   prairieView_Com('zs')
%       tells prairie view to acquire a z stack
%
% METHOD: Given a command or sequence of commands, sends them to prairie
% view over a TCP socket so they will be executed. Returns only when
% prairie view has finished commands. Prairie view monitors port 1236 on
% the broadcasted IP address when opened. It accepts commands following the
% sytax specified in its scripting API. This script keeps the connection to
% the socket open indefintely so it does not need to be reestablished
% (doing so can lead to syncing issues). The draw back of this is if
% prairie view is closed and reopened, then the connection needs to be
% reset by calling prairieView_Com('', "RESET_TCP", 1);
%
% INPUT:
%   cmd_str [str]: a string command or series of commands to be sent to the
%       prairie view computer. Not special delimiters needed, just use
%       spaces as delimiters between arguments, if there is an argument
%       that needs spaces (i.e. a path) use double quotes to encase it.
%   ServerIpAddress[p, str]: ip address of computer running prairie view. If
%       this is null assumes local host.
%   ServerPort [p, double]: port number that is being listened to by the
%       server.
%
% OUTPUT:
%   status [bool]: true if commands are executed, although it is unclear if
%       this will return anything but true, or the proces swill simply hang
%       with an error within prairie.
%
% DEPENDENCIES: nones
%
% AUTHOR: Bruce Corliss, Scientist III
% DATE: 6/16/2012

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

p = inputParser;
p.addRequired('cmd_str', @ischar);
p.addParamValue('ServerIpAddress', [], @ischar);
p.addParamValue('ServerPort', 1236, @ischar);
p.addParamValue('RESET_TCP', 0, @(x) x==1 || x==0);
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

% Reset TCP connection object if specified
if RESET_TCP; setappdata(0, 'micro_tcp_conn', []); end

% Return if input is null
if isempty(cmd_str); status = 1; return; end

% Check to see if microscope communication is disabled
if getappdata(0, 'DISABLE_MICRO_COM'); status=1; return; end

% Establish IP Address
if isempty(ServerIpAddress)
    % Assume host is local, get ip address from system
    host = getLocalIpAddress();
else host = ServerIpAddress;
end

% Port that prairie view is listening to
port = ServerPort;

% Spaces between commands strung together is replaced with character(1)
mult_cmd_str  = cmd_str;
sp_ind = regexp(cmd_str, '\s');
[rs re] = regexp(cmd_str, '".*?["|$]');

% Replace spaces outside of individual arguments with null byte
exarg_space_bw = arrayfun(@(x) any(x>rs' & x<re'), sp_ind);
mult_cmd_str(sp_ind(~exarg_space_bw)) = char(1);
mult_cmd_str = regexprep(mult_cmd_str, [char(1) '*' char(1)],char(1));

% Get base path to code project and current device
proj_path = getappdata(0, 'proj_path');
device_name = getappdata(0, 'device_name');

% Load metadata to determine which port to use
sarr = csvLoad_RowVars([proj_path '/config_clsm/' device_name '.csv']);
port = sarr.micro_com_port;
 
% Assemble socket address.
socket_address = java.net.InetSocketAddress(host,port); 

% Access previous connection if it exists
conn = getappdata(0, 'micro_tcp_conn');
if isempty(conn); socket=[]; else socket = conn.socket; end

%  Attempt to connect if connection does not exist
if isempty(socket) || socket.isClosed() || ~socket.isBound() || ~socket.isConnected()
    % Establish unconnected socket.
    socket = java.net.Socket();
    socket.setSoTimeout(5);
    try socket.connect(socket_address,200);  
    catch ME;  error('%s.m-TCP Connection Failed: %s: %s\n',ME.identifier, ME.message);
    end
    % Get input output streams, put in struct and preserve connection
    conn.socket = socket;
    conn.OutputStream = socket.getOutputStream();
    conn.InputSream = socket.getInputStream();
    setappdata(0, 'micro_tcp_conn', conn);
end
socket.setKeepAlive(true);

% Log command in text file
fid = fopen([proj_path '/tcp_micro_com_log.txt'],'a');
fprintf(fid, '%s TCP[%s:%i]: %s\n', datestr(clock),host, port, mult_cmd_str);
fclose(fid);

%Create output stream, send command to prairie view
outputStream =  conn.OutputStream;
fprintf('>> TCP[%s:%i]: %s', host, port, mult_cmd_str);
outputStream.write(uint8([mult_cmd_str char([13 10])]))
outputStream.flush();

% Create input stream, recieve reply from prairie view
inputStream = conn.InputSream;
% Recieve status message
tcp_end = char([13 10]);
k = 1;
while true
    str = '';
    if inputStream.available() > 0
        for n = 1:inputStream.available()
            str(k+n-1) = inputStream.read();
        end; 
        if strncmp(str(end-1:end), tcp_end, 2); break; end
    end
    pause(.01);
end
% Flush input stream
inputStream.skip(inputStream.available());
%Print to cmd line
str(str==13 | str==10)=[];
fprintf('\t<< %s\n', str);

% Output var
status = strcmp(str, 'DONE');

end

function host = getLocalIpAddress()
% Gets local ip address in winxp or win7 system
strip = @(x) x{1};
[~, ip_str] = system('ipconfig');
if regexp(ip_str, 'IP Address.*:')
    host = strip(regexp(ip_str, 'IP Address.*: ([\d|\.]*)', 'tokens', 'once', 'ignorecase',...
        'dotexceptnewline'));
else % If windows 7 syntax is different
    host = strip(regexp(ip_str, 'IPv4 Address.*: ([\d|\.]*)', 'tokens', 'once', 'ignorecase',...
        'dotexceptnewline'));
end
end


