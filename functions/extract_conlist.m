function cons = extract_conlist(Model)
% Extracts List of Constituents
    [ModName, ~, ~] = rdModFile(Model, 1);
    conList = rd_con(ModName);
    cons = lower(strtrim(cellstr(conList)));
end