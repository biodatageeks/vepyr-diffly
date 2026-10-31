FROM ensemblorg/ensembl-vep:release_116.0
USER root
RUN apt-get update -qq && apt-get install -y -qq libdevel-nytprof-perl > /dev/null && rm -rf /var/lib/apt/lists/*
USER vep
