FROM mambaorg/micromamba:2.0.5
LABEL maintainer="Ricardo R. Pavan (pavan.4@osu.edu)"
LABEL version="1.0.0"
LABEL description="CRESSENT: A comprehensive toolkit for CRESS DNA virus analysis"

# Install CRESSENT from bioconda
RUN micromamba install -y -c conda-forge -c bioconda cressent=1.0.0 && \
    micromamba clean --all --yes

# Set working directory where user data will be mounted
WORKDIR /app

# Set entrypoint to your tool
ENTRYPOINT ["cressent"]

# Default command shows help
CMD ["--help"]
