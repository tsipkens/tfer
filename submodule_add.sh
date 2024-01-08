#!/bin/bash
#
# SUBMODULE_ADD
# A git routine for adding a submodule into a program.
# Author:  Timothy Sipkens, 2020-09-17
#===========================================================#

git submodule add -b main https://github.com/tsipkens/tfer +tfer
git submodule init
