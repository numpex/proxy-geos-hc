# proxy_config_root - the root of the 'configs' folder which includes various platform-dependent .cmake files
export proxy_config_root=`pwd` #'/path/to/proxy-geos-hc' if GUIX_INSTALLED_TPL is ON otherwise /path/to/proxy-geos-hc_tpl
# Name of the config file used to pre-load the cache. The absolute path is ${proxy_config_root}/configs/${config_proxy}
export config_proxy='config_x86_guix.cmake'
#export config_proxy='config_x86_rtx2000.cmake'

# For GUIX_INSTALLED_TPL=OFF, the folder where the TPLs are installed. It is supposed to be a subfolder of $ENV{proxy_config_root}. 
#export install_tpl='installTPLs_O1_C1'
