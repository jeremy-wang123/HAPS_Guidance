function [selected_names, selected_amps, selected_idx] = ...
    select_constituent_amplitudes(glacier_name, Model, mode, specific_constituent)
% SELECT_CONSTITUENT_AMPLITUDES Select tidal constituents by amplitude.
%
% Inputs:
%   glacier_name          Glacier name used by extract_mean_amp
%   Model                 Tidal model path or model identifier
%   mode                  'top', 'top2', 'all', or 'specific'
%   specific_constituent  Constituent name for 'specific' mode
%
% Outputs:
%   selected_names        Cell array of selected constituent names
%   selected_amps         Corresponding mean amplitudes
%   selected_idx          Indices in conlist
%
% Examples:
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'top')
%
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'top2')
%
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'all')
%
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'specific', 'K1')

    if nargin < 4
        specific_constituent = '';
    end

    % Get mean amplitudes
    mean_amps = extract_mean_amp(glacier_name);

    names = fieldnames(mean_amps);
    values = cell2mat(struct2cell(mean_amps));

    % Ensure column vectors
    names = names(:);
    values = values(:);

    % Sort amplitudes from largest to smallest
    [sorted_values, sort_idx] = sort(values,'descend');
    sorted_names = names(sort_idx);

    % Select constituents
    switch lower(mode)

        case 'top'
            selected_names = sorted_names(1);
            selected_amps = sorted_values(1);

        case 'top2'
            n_select = min(2,numel(sorted_names));

            selected_names = sorted_names(1:n_select);
            selected_amps = sorted_values(1:n_select);

        case 'all'
            selected_names = sorted_names;
            selected_amps = sorted_values;

        case 'specific'
            if isempty(specific_constituent)
                error(['A constituent name must be provided ', ...
                       'when mode is ''specific''.']);
            end

            match_idx = find(strcmpi(names,specific_constituent),1);

            if isempty(match_idx)
                error('Constituent "%s" was not found.', ...
                    specific_constituent);
            end

            selected_names = names(match_idx);
            selected_amps = values(match_idx);

        otherwise
            error(['Invalid mode "%s". Choose ''top'', ''top2'', ', ...
                   '''all'', or ''specific''.'],mode);
    end

    % Match selected names to model constituent list
    conlist = extract_conlist(Model);

    % extract_conlist may return a character array
    if ischar(conlist)
        conlist = cellstr(conlist);
    end

    conlist = strtrim(conlist);

    selected_idx = NaN(numel(selected_names),1);

    for k = 1:numel(selected_names)

        idx = find(strcmpi(conlist,selected_names{k}),1);

        if isempty(idx)
            warning('Constituent "%s" was not found in conlist.', ...
                selected_names{k});
        else
            selected_idx(k) = idx;
        end
    end
end