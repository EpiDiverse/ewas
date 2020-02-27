FROM nfcore/base
LABEL authors="Adam Nunn" \
      email="adam.nunn@ecseq.com" \
      description="Docker image containing all requirements for the epidiverse/ewas pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
<<<<<<< HEAD

# you only need to change <template> in this line:
ENV PATH /opt/conda/envs/<template>/bin:$PATH:
=======
ENV PATH /opt/conda/envs/ewas/bin:$PATH:
>>>>>>> refactor for DSL2 changes
