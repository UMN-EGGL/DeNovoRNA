#!/bin/bash

DOCKER_IMAGE="braker2"

function usage() {
	echo "USAGE: $0 PATH"
	echo "	parameters"
	echo "		PATH: path to the directory which"
	echo "		      will become a bind mount under"
	echo "		      /data in the docker container."
	echo "		      All data should be placed here."
}

function run() {
	if [[ $# -ne 1 ]]
	then
		usage
		exit 1
	fi

	# Resolve an absolute path
	BIND_MOUNT_PATH=`realpath "$1"`

	# The punchline
	docker run -it --mount="type=bind,source=$BIND_MOUNT_PATH,destination=/data" $DOCKER_IMAGE
}

run $@
