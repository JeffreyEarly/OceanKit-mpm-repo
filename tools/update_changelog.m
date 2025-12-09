function update_changelog(changelogPath,notes,version)
%% 2) Update CHANGELOG.md (optional notes)
today = datetime('today','Format','yyyy-MM-dd');
entryHeader = sprintf("## [%s] - %s\n", version, char(today));

% Format the NOTES into bullet lines
if strlength(notes) > 0
    noteLines = splitlines(notes);
    bullets = "";
    for k = 1:numel(noteLines)
        line = strtrim(noteLines(k));
        if line ~= ""
            bullets = bullets + "- " + line + newline;
        end
    end
    entryBody = bullets;
else
    entryBody = "- No detailed notes provided.\n";
end

newEntry = entryHeader + entryBody + newline;

defaultHeader = "# Version History" + newline + newline;

if ~isfile(changelogPath)
    % Create new CHANGELOG with header + entry
    changelogText = defaultHeader + newEntry;
else
    oldText = fileread(changelogPath);

    % Detect existing header block (up to the first version heading)
    % Look for the first version heading "## ["
    versionIdx = regexp(oldText, "^## \[", "once", "lineanchors");

    if isempty(versionIdx)
        % No previous version entries → just append after header (or whole file)
        changelogText = oldText;
        if ~endsWith(changelogText, newline)
            changelogText = changelogText + newline;
        end
        changelogText = changelogText + newEntry;

    else
        % Split: header section (everything before first version heading),
        % then the old entries
        headerSection = extractBefore(oldText, versionIdx);
        oldEntries    = extractAfter(oldText, versionIdx - 1);

        % Assemble with newEntry inserted immediately after header
        changelogText = headerSection + newEntry + oldEntries;
    end
end

% Write CHANGELOG.md
fid = fopen(changelogPath, "w");
assert(fid ~= -1, "Could not open CHANGELOG.md for writing");
fwrite(fid, changelogText);
fclose(fid);
end