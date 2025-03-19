#Edit the script according to your setting. The default setting provided with this file assumes that the script is sourced from the proxy-geos-hc folder. 

## proxy_config_root - the root of the 'configs' folder which includes various platform-dependent .cmake files. The default setting would search for a config files (machine specific and config_core.cmake) within proxy-geos-hc 
export proxy_config_root=`pwd` #/path/to/proxy-geos-hc 

## Name of the config file used to pre-load the cache. The absolute path is ${proxy_config_root}/configs/${config_proxy}. In case of failure a trial is given to ${proxy_config_root}/${config_proxy} also
export config_proxy='config_x86_rtx2000.cmake'

# The path to the folder where the TPLs are installed. Not required when installing the TPL from a package manager, i.e BUILD_FROM_TPLMIRROR=OFF
# export install_tpl_dir=`pwd`/../proxy-geos-hc_tpl/installTPLs #/path/to/proxy-geos-hc_tpl/installTPLs
