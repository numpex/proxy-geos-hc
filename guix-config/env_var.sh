# proxy_tpl_dir - the root of the 'configs' folder which includes various platform-dependent .cmake files
export proxy_tpl_dir=`pwd` #'/path/to/proxy-geos-hc'
# Name of the config file used to pre-load the cache. The absolute path is ${proxy_tpl_dir}/configs/${config_tpl}
export config_tpl='config_x86_guix.cmake'

