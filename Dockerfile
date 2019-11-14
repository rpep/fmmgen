FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y python3 python3-dev python3-pip 

RUN python3 -m pip install --upgrade pip
ADD requirements.txt requirements.txt
RUN python3 -m pip install -r requirements.txt
RUN pwd
# install fmmgen package here
ADD . io
RUN python3 -m pip install -e io
WORKDIR /io
CMD ["/bin/bash"]
