FROM ubuntu:20.04
COPY . /src/.

RUN apt-get update
RUN apt-get install -y python3
RUN apt-get install -y python3.8-venv
RUN apt-get install -y bcalm
RUN apt-get install -y jellyfish
RUN apt-get install -y graphviz
RUN apt-get install -y mono-dbg 
RUN apt-get install -y mono-complete
RUN apt-get install -y make
RUN apt-get install -y nuget

# there is a bug in ubuntu 20.04, here, I am solving it using the solution proposed here https://askubuntu.com/questions/1229982/how-to-install-monodevelop-in-20-04-and-get-it-to-build-something
RUN apt-get install -y gnupg ca-certificates
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/ubuntu stable-bionic main" | tee /etc/apt/sources.list.d/mono-official-stable.list
RUN apt-get update
RUN apt-get install -y mono-roslyn

WORKDIR /src
RUN make
