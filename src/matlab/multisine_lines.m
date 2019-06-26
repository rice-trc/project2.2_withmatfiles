function ms = multisine_lines(N, Nt, ms_type)
% Return the excited lines along with harmonics
% Used for time-domain multisine excitation, ie.
%
% phase = 2*pi*rand(N,1);
% fex = @(t) ms.har'*A*cos(2*pi*(1:N)'*f0*t + phase) / sqrt(sum(ms.har));

non_lines = [];
non_odd = [];
non_even = [];

if nargin < 3
    ms_type = 'full';
end
    
% maybe the +1 is wrong. Is is there to convert between DFT lines and lines
% in terms of f0, ie for defining har.
if strcmp(ms_type,'full')
    lines = (1:N)+1;     % % full multisine, eg. excite all lines.
elseif strcmp(ms_type,'odd')
    lines = (1:2:N)+1;    
else  % random odd
    Nblock = 4; % Block length for a random-odd multisine
    lines = (1:2:N); % lines in terms of f0.
    Nblocks = floor(length(lines)/Nblock);
    indices = 1:Nblock:Nblocks*Nblock; % Start indices of blocks ...
    indices = indices+floor(Nblock*rand(1,Nblocks)); % ... plus random number
    lines(indices) = []; % Eliminate selected detection lines  

    % convert to DFT lines! Important!
    lines = lines + 1;
    
    % non-excited lines
    non_lines = 1:Nt/2;
    non_lines(lines) = [];
    non_odd = non_lines(logical(mod(non_lines-1,2)));
    non_even = non_lines(~mod(non_lines-1,2));
end

har = zeros(N,1);
har(lines-1) = 1;

ms.har = har;
ms.lines = lines;
ms.non_lines = non_lines;
ms.non_odd = non_odd;
ms.non_even = non_even;
ms.type = ms_type;

end