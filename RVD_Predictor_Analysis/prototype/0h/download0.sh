#!/bin/bash
# Usage: $ ./download0.sh &> log.download0

for var in syn7197151 syn7203043 syn7205262 syn7209281 syn7209467 syn7209307 syn7202975 syn7209489 syn7203218 syn7205069 syn7185883 syn7203038 syn7203115 syn7443647 syn7208729 syn7208931 syn7207901 syn7188424 syn7187887 syn7186913 syn7209037
do
    echo trying $var...
    synapse get $var
done
