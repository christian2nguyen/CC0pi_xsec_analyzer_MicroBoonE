#!/bin/bash

# Check if exactly two arguments are given
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/dir1 /path/to/dir2"
    exit 1
fi

# Assign arguments to variables
dir1="$1"
dir2="$2"

# Check if the provided arguments are directories
if [ ! -d "$dir1" ]; then
    echo "Error: $dir1 is not a directory."
    exit 1
fi

if [ ! -d "$dir2" ]; then
    echo "Error: $dir2 is not a directory."
    exit 1
fi

# List files in both directories and sort them
files1=$(ls "$dir1" | sort)
files2=$(ls "$dir2" | sort)

# Compare the sorted file lists
if [ "$files1" == "$files2" ]; then
    echo "The directories have the same files."
else
    echo "The directories do not have the same files."
    # Use diff to show differences
    diff <(echo "$files1") <(echo "$files2")
fi

