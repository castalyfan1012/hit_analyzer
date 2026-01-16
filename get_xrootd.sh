#!/bin/bash

while read -r fname; do
    if [ -n "$fname" ]; then
        furl=$(samweb get-file-access-url --schema root "$fname" 2>/dev/null)
        if [ $? -eq 0 ] && [ -n "$furl" ]; then
            echo "$furl"
        else
            echo "Error: Failed to get XRootD URL for $fname" >&2
        fi
    fi
done < filelist_data.txt