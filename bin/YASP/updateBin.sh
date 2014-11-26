#!/bin/sh
#
# Get the latest build result from the 'src/' directory.
#
rsync -avu ../../src/gesis2dDG .
rsync -avu ../../src/gesis3dDG .
rsync -avu ../../src/gesis3dFEM .
rsync -avu ../../src/gesis2dFEM .
