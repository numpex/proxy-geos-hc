#Edit the script according to your setting. The default setting provided with this file assumes that the script is sourced from the proxy-geos-hc folder. 

## The config file which wille be used pre-load the cache / and wrapped by config_proxy-app.cmake (when building the proxyApp)
export config_proxy=`pwd`/configs/'config_x86_rtx2000.cmake'

# The path to the folder where the TPLs are installed. Not required when installing the TPL from a package manager, i.e BUILD_FROM_TPLMIRROR=OFF
# export install_tpl_dir=`pwd`/../proxy-geos-hc_tpl/installTPLs #/path/to/proxy-geos-hc_tpl/installTPLs
