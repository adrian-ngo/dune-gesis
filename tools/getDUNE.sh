#!/bin/sh
#
#
#
echo "Download DUNE modules (tar balls)."
wget http://www.dune-project.org/download/2.3.1/dune-common-2.3.1.tar.gz
wget http://www.dune-project.org/download/2.3.1/dune-geometry-2.3.1.tar.gz
wget http://www.dune-project.org/download/2.3.1/dune-grid-2.3.1.tar.gz
wget http://www.dune-project.org/download/2.3.1/dune-istl-2.3.1.tar.gz
wget http://www.dune-project.org/download/2.3.1/dune-localfunctions-2.3.1.tar.gz
echo "Download DUNE modules done."
#
wget http://dune-project.org/download/pdelab/2.0/dune-typetree-2.3.1.tar.gz
wget http://dune-project.org/download/pdelab/2.0/dune-pdelab-2.0.0.tar.gz
echo "Download DUNE PDELab done."

echo "Start unpacking all."
tar xzvf dune-common-2.3.1.tar.gz
tar xzvf dune-geometry-2.3.1.tar.gz
tar xzvf dune-grid-2.3.1.tar.gz
tar xzvf dune-istl-2.3.1.tar.gz
tar xzvf dune-localfunctions-2.3.1.tar.gz
tar xzvf dune-typetree-2.3.1.tar.gz
tar xzvf dune-pdelab-2.0.0.tar.gz
echo "Unpacking all tar balls done."
