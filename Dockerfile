FROM mambaorg/micromamba:2.0.5
LABEL maintainer="Ricardo R. Pavan (pavan.4@osu.edu)"
LABEL version="1.0.0"

RUN micromamba install -y -n base -c conda-forge -c bioconda genomad==1.0.0 && \
    micromamba clean --all --yes
WORKDIR /app
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "cressent"]