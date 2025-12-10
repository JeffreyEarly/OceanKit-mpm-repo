% mpminstall("/Users/jearly/Documents/ProjectRepositories/class-annotations", Authoring=true); 
% mpminstall("/Users/jearly/Documents/ProjectRepositories/class-docs", Authoring=true); 

pkg = mpmcreate("NetCDF","/Users/jearly/Documents/ProjectRepositories/netcdf");

% Remove the .git folders that were added
paths = string({pkg.Folders.Path});
isGit = startsWith(paths, ".git");
removeFolder(pkg, paths(isGit));
removeFolder(pkg, paths(startsWith(paths, "tools")));
removeFolder(pkg, paths(startsWith(paths, "Documentation")));

pkg.ReleaseCompatibility = ">=R2024b";
pkg.DisplayName = "NetCDF";
pkg.Description = "A class for reading and writing to NetCDF files.";
pkg.Summary =  "A class for reading and writing to NetCDF files.";
pkg.Provider = matlab.mpm.Provider(Name="Jeffrey J. Early", ...
                              Organization="NorthWest Research Associates", ...
                              Email="jearly@nwra.com", ...
                              URL="https://github.com/JeffreyEarly/netcdf");

pkg.addDependency("ClassAnnotations");
pkg.updateDependency("ClassAnnotations","^1.0")