%% setup_oceankit_repos.m
%  1) Make OceanKit-repos folder
%  2) git clone / git pull a list of repos into it
%  3) mpminstall each repo in Authoring mode (editable, in-place)

% ---------- Configuration ----------
baseDir = fullfile("../..", "OceanKit-repos");

% List your repos here
repos = [
    % struct("Name","WaveVortexModelDiagnostics","Url","https://github.com/JeffreyEarly/WaveVortexModelDiagnostics.git","Branch","main")
    struct("Name","class-annotations",          "Url","https://github.com/JeffreyEarly/class-annotations.git",          "Branch","main")
    struct("Name","class-docs",        "Url","https://github.com/JeffreyEarly/class-docs.git",        "Branch","main")
];

% If you use SSH, URLs look like: "git@github.com:JeffreyEarly/OceanKit.git"

% ---------- 1) Ensure base directory exists ----------
if ~isfolder(baseDir)
    mkdir(baseDir);
    fprintf("Created folder: %s\n", baseDir);
else
    fprintf("Folder exists:   %s\n", baseDir);
end

% ---------- 2) Clone or pull ----------
if system("git --version") ~= 0
    error("git does not appear to be available on your PATH. Install git or fix PATH.");
end

for k = 1:numel(repos)
    repoDir = fullfile(baseDir, repos(k).Name);

    if ~isfolder(repoDir)
        % Clone
        fprintf("\n[git] Cloning %s...\n", repos(k).Name);
        cmd = sprintf('git clone --branch "%s" "%s" "%s"', repos(k).Branch, repos(k).Url, repoDir);
        [status, out] = system(cmd);
        if status ~= 0
            error("git clone failed for %s:\n%s", repos(k).Name, out);
        end
    else
        % Pull (only if it looks like a git repo)
        if isfolder(fullfile(repoDir, ".git"))
            fprintf("\n[git] Updating %s...\n", repos(k).Name);

            % Fetch + checkout branch + pull (more robust than plain pull)
            cmd = sprintf('git -C "%s" fetch --all --prune', repoDir);
            [status, out] = system(cmd);
            if status ~= 0, error("git fetch failed for %s:\n%s", repos(k).Name, out); end

            cmd = sprintf('git -C "%s" checkout "%s"', repoDir, repos(k).Branch);
            [status, out] = system(cmd);
            if status ~= 0, error("git checkout failed for %s:\n%s", repos(k).Name, out); end

            cmd = sprintf('git -C "%s" pull', repoDir);
            [status, out] = system(cmd);
            if status ~= 0, error("git pull failed for %s:\n%s", repos(k).Name, out); end
        else
            warning("%s exists but is not a git repo (no .git folder). Skipping update.", repoDir);
        end
    end

    % ---------- 3) Install in authoring mode ----------
    fprintf("[mpm] Installing (Authoring=true): %s\n", repoDir);

    pkgFile = fullfile(repoDir, "mpackage.json");
    if ~isfile(pkgFile)
        error("Expected %s to contain mpackage.json, but it was not found: %s", repos(k).Name, pkgFile);
    end

    try
        pkg = mpminstall(repoDir, Authoring=true);
        fprintf("[mpm] Done: %s (installed/updated)\n", pkg.Name);
    catch ME
        fprintf(2, "[mpm] FAILED for %s\n", repos(k).Name);
        rethrow(ME);
    end
end

fprintf("\nAll repos processed.\n");