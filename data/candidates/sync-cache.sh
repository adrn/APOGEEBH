#!/bin/bash

rsync -zvr --max-size=1000m --exclude "*~" --exclude "*.npy" \
perseus:/tigress/adrianp/projects/apogeebh/data/candidates/cache ~/projects/apogeebh/data/candidates
