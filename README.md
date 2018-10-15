# rnaseq-pipeline

## Prerequisites

1. Install [docker CE](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

## Building

```
docker image build -t rnaseq-pipe .
```

## Running

## Development

```
docker run -it --rm -v $PWD:/home/ rnaseq-pipe
```

This maps the current directory to a folder inside the container. This makes it
so that there is no need to re-build the image upon making source code changes.
