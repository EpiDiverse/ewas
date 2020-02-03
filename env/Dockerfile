FROM nfcore/base
LABEL authors="Firstname Surname" \
      email="firstname.surnname@domain.com" \
      description="Docker image containing all requirements for the epidiverse/template pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# you only need to change <template> in this line:
ENV PATH /opt/conda/envs/<template>/bin:$PATH:
