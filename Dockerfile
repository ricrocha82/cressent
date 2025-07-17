cat > Dockerfile << 'EOF'

FROM mambaorg/micromamba:2.0.5
LABEL maintainer="Ricardo R. Pavan (pavan.4@osu.edu)"
LABEL version="1.0.0"
LABEL description="CRESSENT: A comprehensive toolkit for CRESS DNA virus analysis"

# Install mamba for faster package management
RUN conda install -c conda-forge mamba -y

# Install from bioconda (recommended)
RUN mamba install -y -c conda-forge -c bioconda cressent=1.0.0

# Clean up to reduce image size
RUN mamba clean --all --yes && \
    rm -rf /var/lib/apt/lists/*

# Set working directory where user data will be mounted
WORKDIR /app

# Set entrypoint to your tool
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "cressent"]

# Default command shows help
CMD ["--help"]

EOF
