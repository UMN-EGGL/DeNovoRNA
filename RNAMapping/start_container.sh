#!/bin/bash
docker run \
    --mount type=bind,source=/project,target=/project \
    --mount type=bind,source=/scratch,target=/scratch \
    --mount type=bind,source=/home/rob/Codes/DeNovoRNA,target=/root/Codes/DeNovoRNA \
    --volume star_muscle_fat:/output \
    -it rnamap
