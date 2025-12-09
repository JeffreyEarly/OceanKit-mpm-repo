function ci_release(options)
arguments
    options.rootDir = ".."
    options.bumpType = "none"
    options.notes string = ""
    options.shouldBuildWebsiteDocumentation = false
    options.shouldPackageForDistribution = false
    options.dist_folder = "dist"
    options.excluded_dist_folders = [".git", ".github", "docs", "tools", "Documentation"]
end
%CI_RELEASE CI entry point for MPM release.
%   CI_RELEASE(options) where
%       options.bumpType is "patch", "minor", or "major". If it is left
%
%   Steps:
%     1) Bump version in resources/mpackage.json
%     2) Run custom documentation build
%     3) Export package root to dist/<name>-<version> for MPM repo

if ~isfield(options,"rootDir")
    options.rootDir = pwd;
end

%% bump the version number, if requested
mpkgPath = fullfile(options.rootDir, "resources", "mpackage.json");
if ~isfile(mpkgPath)
    error("ci_release:mpackageNotFound", ...
        "Could not find mpackage.json at %s", mpkgPath);
end

% 1) Read and bump version (semantic x.y.z)
txt = fileread(mpkgPath);
cfg = jsondecode(txt);

if ~isfield(cfg,"version")
    error("ci_release:noVersionField", ...
        "mpackage.json does not have a 'version' field.");
end
if ~isfield(cfg,"name")
    error("ci_release:noNameField", ...
        "mpackage.json does not have a 'name' field.");
end

oldVer = string(cfg.version);
tokens = regexp(oldVer, "^(\d+)\.(\d+)\.(\d+)$", "tokens", "once");
if isempty(tokens)
    error("ci_release:badVersion", ...
        "Existing version '%s' is not of the form x.y.z", oldVer);
end
major = str2double(tokens{1});
minor = str2double(tokens{2});
patch = str2double(tokens{3});

if options.bumpType ~= "none"
    switch options.bumpType
        case "major"
            major = major + 1;
            minor = 0;
            patch = 0;
        case "minor"
            minor = minor + 1;
            patch = 0;
        case "patch"
            patch = patch + 1;
    end

    newVer = sprintf('%d.%d.%d', major, minor, patch);
    cfg.version = newVer;

    fprintf('Bumping version: %s -> %s (%s)\n', oldVer, newVer, options.bumpType);

    % Write back pretty JSON (R2021b+ supports PrettyPrint)
    jsonStr = jsonencode(cfg, PrettyPrint=true);
    fid = fopen(mpkgPath, 'w');
    assert(fid ~= -1, "Could not open mpackage.json for writing");
    fwrite(fid, jsonStr);
    fclose(fid);

    %% If we bumped the version, and were handed notes, record that
    changelogPath = fullfile(options.rootDir, "CHANGELOG.md");
    update_changelog(changelogPath,options.notes,newVer);
else
    newVer = sprintf('%d.%d.%d', major, minor, patch);
end

%% 2) Run your custom documentation build
if options.shouldBuildWebsiteDocumentation
    % Replace this with your actual doc build entry point
    % e.g. waveVortexDiagnostics_build_docs, or build_docs
    if exist("build_website_documentation","file")
        fprintf('Running documentation builder\n');
        build_website_documentation(rootDir=options.rootDir);
    end
end

%% 3) Export package root to dist/<name>-<version> for MPM repo

if options.shouldPackageForDistribution == true
    pkgName = string(cfg.name);
    distDir = fullfile(options.rootDir, options.dist_folder);
    if ~isfolder(distDir)
        mkdir(distDir);
    end

    pkgFolderName = pkgName + "-" + newVer;
    targetRoot    = fullfile(distDir, pkgFolderName);

    % Clean any stale output
    if isfolder(targetRoot)
        rmdir(targetRoot, "s");
    end

    fprintf('Exporting package root to %s\n', targetRoot);
    copyfile(options.rootDir, targetRoot);

    % Strip CI-only junk from the exported package
    % (best-effort: ignore errors if these don't exist)
    for iFolder = 1:length(options.excluded_dist_folders)
        try
            rmdir(fullfile(targetRoot, options.excluded_dist_folders(iFolder)), "s");
        catch
        end
    end

    try
        rmdir(fullfile(targetRoot, "dist"), "s");
    catch
    end

    % Write a small metadata file for the GitHub Action
    metaPath = fullfile(distDir, "mpm_release_metadata.txt");
    fid = fopen(metaPath, "w");
    assert(fid ~= -1, "Could not open metadata file for writing");
    fprintf(fid, "NAME=%s\nVERSION=%s\nFOLDER=%s\n", ...
        pkgName, newVer, pkgFolderName);
    fclose(fid);

    fprintf('ci_release complete: %s %s\n', pkgName, newVer);
end
end