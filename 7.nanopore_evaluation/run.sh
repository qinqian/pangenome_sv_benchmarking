#!/bin/bash 

load_env() {
    eval "$(/homes6/alvin/miniforge3/bin/conda shell.bash hook)"
}

main() {
    load_env
}

main
