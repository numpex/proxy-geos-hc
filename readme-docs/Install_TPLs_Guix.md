# Getting Guix package manager 

## Guix installation 
Install Guix  following the steps described in the [Guix documentation](https://guix.gnu.org/manual/fr/guix.fr.html#Installation). We recommand to install Guix from the [script installer](https://guix.gnu.org/manual/en/html_node/Binary-Installation.html#index-installer-script) in order to have the latest updates of the `guix-daeomon`. 

## Getting started - guix shell 
There are two options to build an isolated development environment. One involves using `guix shell --container`, the other `guix shell --pure`. In the following, we describe the latest approach. The `guix shell` command allows to create a child environment within a parent environment. Getting back to the parent environment can be achieved running `exit` or `Ctrl+D`. The `--pure` option is useful to clear the environment variables of the parent environment. Therefore enhancing the isolation.  
The command 
```
guix shell --pure <package_list>
```
 will create a development environment with the package list  `<package_list>`. The packages to be installed can be listed in a (Scheme) file, so called manifest file, which can be passed to the `guix shell` command with the option `-m`:
```
guix shell --pure -m manifest.scm
```
A manifest file can be created passing the option `--export-manifest` to guix shell command and redirecting the output in a file
```
guix shell --pure <package_list> --export-manifest > manifest.scm
```
# Creating the  development environment
## Development environment building the TPLs and the proxyApp's main 
The manifest file `proxy-geos-hc/guix-manifest/manifest_proxy_prerequisites.scm` is useful to create a development environment for both the TPLs and the proxyApp.
## Development environment building the proxyApp's main
To get a development environment with both the prerequisites and the TPLs packages, several flavours of manifest file are provided in the folder `guix-manifest`. 
The TPLs required to build the proxyApp depends on the abstraction library to be used (`Raja` or `Kokkos`). For `raja`, the `camp`, `chai` and `umpire` TPLs are also required with their specific configuration depending on the backend programming models enabled.  
### guix search `<tpl_name>`
Use `guix search` command to find out if a package is available within the guix channels. The command reads   
``` 
guix search <tpl_name>* | grep "name: <tpl_name>" 
```
`<tpl_name>` being the placeholder for one of the required TPLs (`raja, kokkos,camp, chai, umpire`). The `*` is added to find all the package names starting with the string  <tpl_name>. As an illustration, for `kokkos`, one would run the command 
```
guix search kokkos* | grep "name: kokkos"
```

### Local channel within the `guix-config` folder
The available packages may not suit the configuration you are considering. With `guix`, there is flexibility to change the arguments for the `cmake build-system` of the library (`--with-configure-flag` option), specifying a different device `DCMake_CUDA_ARCHITECTURES` or host architecture (`--tune` option). It is also possible to `build` for a specific `commit` version (`--with-commit` option) or disable the tests associated to the packages build (`--without-tests` option). The `proxy-geos-hc` repository provides with an internal channel with several package variants. You may consider them using the option `-L` in the guix commands. Therefore, one could run instead
```
guix search -L proxy-geos-hc/guix-config <tpl_name> | grep "name: <tpl_name>"
```

### Create your manifest file 
Create a new manifest file which accounts for the specific TPLs that are required
and the prerequisites
```
guix shell -L proxy-geos-hc/guix-config --pure -m proxy-geos-hc/guix-manifest/manifest_proxy_prerequisites.scm <tpls_packages_list> --export-manifest > geos-hc/guix-manifest/manifest_tpls_<string>_<string>.scm
```
<details><summary>A tip for creating a manifest with package transformations</summary>    
One of the simplest way for creating a manifest with some package transformations for different packages could be to create a manifest file for each of the packages and then gather them in a unique manifest file. Indeed one could run the `guix shell` command passing the `-m` option a multiple time.

```
guix shell -m manifest_1.scm -m manifest_2.scm -m ... --export-manifest manifest_merged.scm
```
It worth to notice that the `--pure` option not required there, since the aim with the previous commandline is simply to create the manifest file and not to create the development environment. 
</details> 
### Building the development environment for the proxyApp's main
Now build the development environment by running the commandline 
```
guix shell -L proxy-geos-hc/guix-config --pure -m proxy-geos-hc/guix-manifest/manifest_tpls_<raja or kokkos>_<string>.scm
```
There are several examples of manifest files provided in the folder guix-manifest.  

Once the development environment has been created, continue with the [build and install of the proxyApp.](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/tree/dev_docs/guix?ref_type=heads#step-3-build-and-install-the-proxyapp)

<!-- # Build the proxyApp's main 
To build the proxyApp within the development environment, run the following commandline 
```
cmake  -DBUILD_FROM_TPLMIRROR=OFF -DCMAKE_BUILD_TYPE=RELEASE <KOKKOS_RAJA_OMP> -C configs/config_proxy-app.cmake -B build -DCMAKE_INSTALL_PREFIX=install -S .
```
Setting `-DBUILD_FROM_TPLMIRROR=OFF` allows the CMake build system to consider finding the installed packages without any path information from the user.    -->