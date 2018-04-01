FROM continuumio/miniconda3:latest

# Update conda
RUN conda update -n base conda --yes
RUN conda update --all --yes

# Download poreduck
RUN git clone -b yield_plot https://github.com/alexiswl/poreduck.git
WORKDIR ./poreduck

# Install required packages
RUN conda install --file requirements.txt --yes

# Re-update conda
RUN conda update --all --yes

# Install poreduck using pip
RUN pip install -e .

# Change to /data directory
WORKDIR /data

# Copy the entry point for the user
COPY ./docker-entrypoint.sh /

# Change user
USER docker

# Set the entrypoint to be 'poreduck'
ENTRYPOINT ["/docker-entrypoint.sh"]
