# Builder stage:
FROM ubuntu:20.04
RUN apt-get update && apt-get install -y --fix-missing python3 && apt-get install -y python3.8-venv \
	&& apt-get install -y bcalm \
	&& apt-get install -y jellyfish \
	&& apt-get install -y graphviz \
	&& apt-get install -y mono-dbg \
	&& apt-get install -y mono-complete \
	&& apt-get install -y make \
	&& apt-get install -y nuget
	
# Create the virtual environment
RUN python3 -m venv /src/venv
ENV PATH=/src/venv/bin:$PATH

# Install Python dependencies
WORKDIR /src
COPY PYTHON_REQUIREMENTS.txt .
RUN . /src/venv/bin/activate && pip3 install wheel && pip3 install -r PYTHON_REQUIREMENTS.txt


# there is a bug in ubuntu 20.04, here, I am solving it using the solution proposed here https://askubuntu.com/questions/1229982/how-to-install-monodevelop-in-20-04-and-get-it-to-build-something
RUN apt-get install -y gnupg ca-certificates \ 
	&& apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF \
	&& echo "deb https://download.mono-project.com/repo/ubuntu stable-bionic main" | tee /etc/apt/sources.list.d/mono-official-stable.list \ 
	&& apt-get update \
	&& apt-get install -y mono-roslyn

# Final stage:
COPY ./alignment /src/alignment/.
COPY ./decomposition /src/decomposition/.
COPY ./example_input /src/example_input/.
COPY ./peak_calling /src/peak_calling/.
COPY ./pseudoassembly /src/pseudoassembly/.
COPY ./makefile /src/.
COPY ./protopaperfly /src/.
COPY ./README.md /src/.

RUN cd /src/ \ 
	&& make build && make construct
ENV PATH /src/:$PATH
WORKDIR /src/

