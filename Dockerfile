FROM ubuntu:22.04

RUN apt update -qq \
    && apt -y --no-install-recommends install \
        build-essential \
        wget

# Python
RUN apt -y --no-install-recommends install \
        python3-dev \
        python3-pip

RUN pip3 install pandas biopython snakemake \
        -i https://pypi.tuna.tsinghua.edu.cn/simple \
        --default-timeout=100

# R
RUN apt -y --no-install-recommends install \
        software-properties-common \
        dirmngr \
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC \
        apt -y --no-install-recommends install tzdata \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
        tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt -y --no-install-recommends install r-base

# Tidyverse
RUN apt -y --no-install-recommends install \
        libxml2-dev \
        libssl-dev \
        libfontconfig1-dev \
        libcurl4-openssl-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
    && Rscript -e 'install.packages("tidyverse")'
        
# Bioconductor
RUN Rscript -e 'install.packages("BiocManager")' \
    && Rscript -e 'BiocManager::install(version = "3.18")'

# Biostrings
RUN apt -y --no-install-recommends install zlib1g-dev libcurl4-gnutls-dev \
    && Rscript -e 'BiocManager::install("Biostrings", update = TRUE, ask = FALSE)'

# fastp
RUN wget -O fastp http://opengene.org/fastp/fastp.0.23.4 \
    && chmod +x fastp \
    && mv fastp /usr/local/bin/

# samtools
RUN apt -y --no-install-recommends install samtools

WORKDIR /isomir
COPY core ./core
COPY script ./script
COPY Snakefile ./

CMD ["snakemake", "--core", "4", "--forceall", "--quiet", "all"]