#!/bin/sh
if [ $(uname) == "Darwin" ]; then
    open video/output.mp4
else
    xdg-open video/output.mp4
fi
