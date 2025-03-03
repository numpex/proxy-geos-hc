# proxy_config_root - the root of the 'configs' folder which includes various platform-dependent .cmake files
export proxy_config_root=`pwd` #'/path/to/proxy-geos-hc' if BUILD_FROM_TPLMIRROR is OFF otherwise /path/to/proxy-geos-hc_tpl
# Name of the config file used to pre-load the cache. The absolute path is ${proxy_config_root}/configs/${config_proxy}
export config_proxy='config_x86_guix.cmake'
#export config_proxy='config_x86_rtx2000.cmake'

# For BUILD_FROM_TPLMIRROR=ON, the build and install folders for the TPLs. It is supposed to be a subfolder of $ENV{proxy_config_root}. 
#export build_tpl='buildTPLs_O1_C1'
#export install_tpl='installTPLs_O1_C1'
