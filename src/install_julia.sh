#!/bin/bash
if [ -d "/shared" ]
then
    cd /shared
else
    cd
fi
wget https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.1-linux-x86_64.tar.gz
tar xvfa julia-1.4.1-linux-x86_64.tar.gz
rm julia-1.4.1-linux-x86_64.tar.gz
if [ -d "/shared" ]
then
    echo PATH=\$PATH:/shared/julia-1.4.1/bin/ >> ~/.bashrc
else
    echo PATH=\$PATH:~/julia-1.4.1/bin/ >> ~/.bashrc
fi
source ~/.bashrc
julia ~/Clements-Denolle-2020/src/add-packages.jl
cd
