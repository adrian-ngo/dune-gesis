#!/bin/bash
#
#
# Download all DUNE modules necessary for DUNE-GESIS via git:
#
#
git clone -b releases/2.3 http://git.dune-project.org/repositories/dune-common
git clone -b releases/2.3 http://git.dune-project.org/repositories/dune-geometry
git clone -b releases/2.3 http://git.dune-project.org/repositories/dune-grid
git clone -b releases/2.3 http://git.dune-project.org/repositories/dune-istl
git clone -b releases/2.3 http://git.dune-project.org/repositories/dune-localfunctions
git clone -b releases/2.3 http://git.dune-project.org/repositories/dune-typetree
git clone -b releases/2.0 http://git.dune-project.org/repositories/dune-pdelab
#
#
# Download the DUNE module for ALUGrid (optional for DUNE-GESIS) via git:
#
# git clone -b releases/2.3 http://users.dune-project.org/repositories/projects/dune-alugrid.git
