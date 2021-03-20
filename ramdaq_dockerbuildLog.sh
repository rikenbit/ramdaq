#!/bin/bash

echo 'DOCKER_BUILD MODE =' "${1:-FALSE}"
echo 'DOCKER_PUSH MODE =' "${2:-FALSE}"
echo 'GIT_COMMIT MODE =' "${3:-FALSE}"

DOCKER_BUILD="${1:-FALSE}"
DOCKER_PUSH="${2:-FALSE}"
GIT_COMMIT="${3:-FALSE}"

if [ "$DOCKER_BUILD" = "TRUE" ];
then
    echo ">>> The docker build will be executed..."
fi

if [ "$DOCKER_PUSH" = "TRUE" ];
then
    echo ">>> The docker image will be automatically pushed..."
fi

if [ "$GIT_COMMIT" = "TRUE" ];
then
    echo ">>> The ramdaq source code will be automatically committed..."
fi

IMAGE_NAME=myoshimura080822/ramdaq
IMAGE_TAG=v1.1

DATE=$(date '+%Y%m%d-%H%M%S')
REPO_COMMIT=`git log -n 1 --pretty=format:%H | cut -c1-7`
YAML_COMMIT=`git log -n 1 --pretty=format:%H -- environment.yml | cut -c1-7`

# Make log directory
mkdir -p log

# Define log file name
LOGFILE="log/log_docker_${IMAGE_TAG}_${DATE}.txt"

# Record date
echo "DATE=$DATE" > ${LOGFILE}

# Record Docker image tag
echo "IMAGE_TAG=${IMAGE_TAG}" >> ${LOGFILE}

# Record the current commit ID
echo "REPO_COMMIT=${REPO_COMMIT}" >> ${LOGFILE}
echo "YAML_COMMIT=${YAML_COMMIT}" >> ${LOGFILE}

if [ "$DOCKER_BUILD" = "TRUE" ];
then
    # Build Docker image
    cmd="docker build -t ${IMAGE_NAME}:${YAML_COMMIT} ."
    echo -e "\nDocker log:" >> ${LOGFILE}
    echo $cmd >> ${LOGFILE}
    $cmd
    
    # Tag Docker image
    cmd="docker tag ${IMAGE_NAME}:${YAML_COMMIT} ${IMAGE_NAME}:${IMAGE_TAG}"
    echo $cmd >> ${LOGFILE}
    $cmd
fi

if [ "$DOCKER_PUSH" = "TRUE" ];
then
    # Push Docker image
    cmd="docker push ${IMAGE_NAME}:${IMAGE_TAG}"
    echo $cmd >> ${LOGFILE}
    $cmd
fi

if [ "$GIT_COMMIT" = "TRUE" ];
then
    # git add and commit
    git add .
    git commit -m "build docker img ${YAML_COMMIT} as ${IMAGE_TAG}"
fi
