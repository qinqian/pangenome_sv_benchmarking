#!/bin/bash
squeue -u $USER | grep 7 | awk '{print $1}' | xargs -n 1 scancel
