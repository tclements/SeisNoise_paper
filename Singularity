Bootstrap: library
From: ubuntu:16.04
Stage: build

%files
    ./src/add-packages.jl add-packages.jl

%environment
    export PATH=/julia-1.1.1/bin:$PATH
    export LD_LIBRARY_PATH=/julia-1.1.1/lib:/julia-1.1.1/lib/julia:$LD_LIBRARY_PATH
    export LC_ALL=C

%post
    apt-get update
    apt-get -y install wget tar
    apt-get -y install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6

    # install julia
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.1-linux-x86_64.tar.gz
    tar  -vzxf julia-1.1.1-linux-x86_64.tar.gz
    rm -f julia-1.1.1-linux-x86_64.tar.gz
    ln -s /julia-1.1.1/bin/julia /usr/local/bin/julia

    # install julia packages 
    julia add-packages.jl
    
%runscript
    exec "$@"

%test
    grep -q NAME=\"Ubuntu\" /etc/os-release
    if [ $? -eq 0 ]; then
        echo "Container base is Ubuntu as expected."
    else
        echo "Container base is not Ubuntu."
    fi

%labels
    Author thclements@g.harvard.edu
    Version v0.0.1

%help
    This is a demo container for running Clements and Denolle, 2019 on AWS. 