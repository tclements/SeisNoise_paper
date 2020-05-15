#!/bin/bash
bash -ci "$(curl -fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
source ~/.bashrc
julia ~/SeisNoise_paper/src/add-packages.jl
