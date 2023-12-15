# Indicate the Gurobi reference image
FROM gurobi/python:11.0.0_3.11 

# Set the application directory
WORKDIR /app

# Update the OS and install dev tools
## for apt to be noninteractive
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt update -qqy
RUN apt-get install -qqy apt-utils

RUN echo 'tzdata tzdata/Areas select Etc' | debconf-set-selections; \
    echo 'tzdata tzdata/Zones/Etc select UTC' | debconf-set-selections; \
    apt-get install -qqy --no-install-recommends tzdata

## preesed tzdata, update package index, upgrade packages and install needed software
RUN apt-get install -qqy \
    build-essential \
    ca-certificates \
    cmake \
    dirmngr \
    make \
    util-linux \
    wget 
RUN apt-get clean

# Install the application and its dependencies
ADD requirements.txt /app/
RUN pip install -r requirements.txt

# Prepare the data volumes
RUN mkdir /data /results
VOLUME ["/data"]
VOLUME ["/results"]
VOLUME ["/config"]

# Install the app only at the end, so dockerfile development is faster
ENV PATH=$PATH:/app
ADD . /app
# Command used to start the application
# CMD ["/app/pangeblock"]
CMD ["/usr/local/bin/snakemake", "-s", "/app/pangeblocks.smk", "-c32", "--configfile=/config/params.yml"]
# CMD ["touch /results/results.txt"]
