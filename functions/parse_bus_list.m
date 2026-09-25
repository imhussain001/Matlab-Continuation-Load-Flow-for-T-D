function buses = parse_bus_list(str, busTypes)
%PARSE_BUS_LIST  Parse a user bus selection such as '14', '4 9 14', '10-14, 18', 'all', 'pq'.
%
%   buses = parse_bus_list(str, busTypes)
%
%   busTypes : N x 1 bus types (3 slack, 2 PV, 0 PQ); N defines valid buses
%   'all'    : every bus except the slack
%   'pq'     : load (PQ) buses only
%   'a-b'    : inclusive range
%   Separators: spaces, commas or semicolons. Duplicates are removed and
%   the original order kept. Throws parse_bus_list:Invalid with a readable
%   message for anything else.

N = numel(busTypes);
str = lower(strtrim(string(str)));
if str == "all"
    buses = find(busTypes ~= 3).';
    return
elseif str == "pq"
    buses = find(busTypes == 0).';
    return
end

tokens = split(regexprep(str, '\s*-\s*', '-'), {' ', ',', ';'});
tokens = tokens(tokens ~= "");
if isempty(tokens)
    error('parse_bus_list:Invalid', 'No bus numbers given.');
end
buses = zeros(1,0);
for tok = tokens.'
    parts = split(tok, '-');
    nums = str2double(parts);
    if any(isnan(nums)) || any(nums ~= round(nums)) || numel(nums) > 2
        error('parse_bus_list:Invalid', '"%s" is not a bus number or range (use e.g. 14 or 10-14).', tok);
    end
    if numel(nums) == 2
        nums = nums(1):sign(nums(2)-nums(1)+eps):nums(2);
    end
    buses = [buses, nums(:).']; %#ok<AGROW>
end
bad = buses(buses < 1 | buses > N);
if ~isempty(bad)
    error('parse_bus_list:Invalid', 'Bus %s does not exist (valid buses: 1-%d).', mat2str(bad), N);
end
buses = unique(buses, 'stable');
end
