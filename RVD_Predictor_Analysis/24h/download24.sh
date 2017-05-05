#!/bin/bash
# Usage: $ ./download24.sh &> log.download24

for var in syn7202904 syn7203292 syn7205265 syn7209284 syn7209466 syn7209400 syn7208033 syn7209490 syn7203219 syn7205070 syn7207859 syn7203145 syn7203150 syn7443645 syn7208730 syn7208936 syn7208928 syn7207898 syn7188424 syn7187888 syn7186914 syn7209039                  
do
    echo trying $var...
    synapse get $var
done
