% mpminstall("/Users/jearly/Documents/ProjectRepositories/class-annotations", Authoring=true); 
% mpminstall("/Users/jearly/Documents/ProjectRepositories/class-docs", Authoring=true); 

pkg = mpmcreate("AlongTrackSimulator","/Users/jearly/Documents/ProjectRepositories/AlongTrackSimulator");

% Remove the .git folders that were added
paths = string({pkg.Folders.Path});
isGit = startsWith(paths, ".git");
removeFolder(pkg, paths(isGit));
removeFolder(pkg, paths(startsWith(paths, "tools")));
removeFolder(pkg, paths(startsWith(paths, "Examples")));

pkg.ReleaseCompatibility = ">=R2024b";
pkg.DisplayName = "AlongTrackSimulator";
pkg.Description = "Ground track simulator for altimetry missions.";
pkg.Summary =  "Ground track simulator for altimetry missions.";
pkg.Provider = matlab.mpm.Provider(Name="Jeffrey J. Early", ...
                              Organization="NorthWest Research Associates", ...
                              Email="jearly@nwra.com", ...
                              URL="https://github.com/satmapkit/AlongTrackSimulator");

% pkg.addDependency("ClassAnnotations");
% pkg.updateDependency("ClassAnnotations","^1.0")