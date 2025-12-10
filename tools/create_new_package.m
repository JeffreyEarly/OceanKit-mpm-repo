% mpminstall("/Users/jearly/Documents/ProjectRepositories/class-annotations", Authoring=true); 
% mpminstall("/Users/jearly/Documents/ProjectRepositories/class-docs", Authoring=true); 

pkg = mpmcreate("WaveVortexModelDiagnostics","/Users/jearly/Documents/ProjectRepositories/wave-vortex-model-diagnostics");

% Remove the .git folders that were added
paths = string({pkg.Folders.Path});
isGit = startsWith(paths, ".git");
removeFolder(pkg, paths(isGit));

pkg.ReleaseCompatibility = ">=R2024b";
pkg.DisplayName = "WaveVortexModel Diagnostics";
pkg.Description = "A class for quickly analyzing output from the wave-vortex model, with a focus on examining the energy fluxes from forcing and triads.";
pkg.Summary = "Analyze output from the WaveVortexModel";
pkg.Provider = matlab.mpm.Provider(Name="Jeffrey J. Early", ...
                              Organization="NorthWest Research Associates", ...
                              Email="jearly@nwra.com", ...
                              URL="https://github.com/JeffreyEarly/distributions");

pkg.addDependency("SplineCore");
pkg.updateDependency("SplineCore","^1.3")