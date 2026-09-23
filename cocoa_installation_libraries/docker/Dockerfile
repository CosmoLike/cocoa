FROM ubuntu:22.04

ARG NB_USER="whovian"
ARG NB_UID=1000
ARG NB_GID=1000
# Cocoa version whose conda-lock yml file seeds the cocoa environment
ARG COCOA_TAG="v4.11.8"

# Fix DL4006
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install -yq --fix-missing software-properties-common  \
    apt-utils build-essential \
  && apt-get update \
  && apt-get install -yq --fix-missing --no-install-recommends \
    git autoconf dpkg xz-utils pkg-config \
  && apt-get install -yq --fix-missing --no-install-recommends libgflags-dev \
    nano make wget curl libgtk2.0-dev libglib2.0-0 libfreetype6-dev libleveldb-dev \
    tzdata unzip grep libsnappy-dev sed sudo locales fonts-dejavu \
    libtool automake zlib1g-dev libcurl4-gnutls-dev libssl-dev run-one emacs sed vim \
    texlive texlive-xetex  texlive-latex-extra texlive-plain-generic \
    texlive-fonts-recommended dvipng cm-super ghostscript dvipng less \
# ------------------------------------------------------------------------
  && ln -s /bin/tar /bin/gtar \
  && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen \
  && locale-gen \
# ------------------------------------------------------------------------
  && curl -L "https://github.com/krallin/tini/releases/download/v0.19.0/tini_0.19.0.deb" > tini.deb \
  && dpkg -i tini.deb \
  && rm tini.deb \
# --------------------------------------------------------------------------
# Install OpenSSH for MPI to communicate between containers
  && apt-get update \
  && apt-get install -y --no-install-recommends openssh-client openssh-server \
  && mkdir -p /var/run/sshd \
# --------------------------------------------------------------------------
# Allow OpenSSH to talk to containers without asking for confirmation
  && cat /etc/ssh/ssh_config | grep -v StrictHostKeyChecking > /etc/ssh/ssh_config.new \
  && echo "    StrictHostKeyChecking no" >> /etc/ssh/ssh_config.new \
  && mv /etc/ssh/ssh_config.new /etc/ssh/ssh_config \
# --------------------------------------------------------------------------
  && apt-get clean \
  && apt-get autoremove --yes \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/doc/ \
  && rm -rf /usr/share/man/

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# GLOBAL ENV VARIABLES

ENV CONDA_DIR=/home/$NB_USER/conda \
    SHELL=/bin/bash \
    NB_USER=$NB_USER \
    NB_UID=$NB_UID \
    NB_GID=$NB_GID \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    HOME=/home/$NB_USER \
    MAKE_NUM_THREADS=4 \
    OMP_PROC_BIND=close \
    OMP_NUM_THREADS=3 \
    GOSU_VERSION=1.16

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Copy a script that we will use to correct permissions
# after running certain commands
COPY fix-permissions /usr/local/bin/fix-permissions
RUN chmod a+rx /usr/local/bin/fix-permissions \
  && sed -i 's/^#force_color_prompt=yes/force_color_prompt=yes/' /etc/skel/.bashrc
# --------------------------------------------------------------------------
# https://github.com/jupyter/docker-stacks/pull/826
# Enable prompt color in the skeleton .bashrc before creating the NB_USER

# --------------------------------------------------------------------------
RUN curl -OL "https://github.com/google/protobuf/releases/download/v3.0.0/protoc-3.0.0-linux-x86_64.zip" \
  && unzip protoc-3.0.0-linux-x86_64.zip -d proto3 \
  && mv proto3/bin/* /usr/local/bin \
  && mv proto3/include/* /usr/local/include \
  && rm -rf proto3 protoc-3.0.0-linux-x86_64.zip
# --------------------------------------------------------------------------
# INSTALL GOSU
RUN dpkgArch="$(dpkg --print-architecture | awk -F- '{ print $NF }')" \
  && curl -sSL "https://github.com/tianon/gosu/releases/download/${GOSU_VERSION}/gosu-${dpkgArch}" > /usr/bin/gosu \
  && chmod 755 /usr/bin/gosu \
  && gosu nobody true \
# --------------------------------------------------------------------------
# INSTALL wait-for-it.sh
  && curl -sSL https://github.com/vishnubob/wait-for-it/raw/master/wait-for-it.sh > /usr/bin/wait-for-it.sh \
  && chmod 755 /usr/bin/wait-for-it.sh \
# --------------------------------------------------------------------------
  && ldconfig

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Create NB_USER wtih name jovyan user with UID=1000 and in the 'users' group
# and make sure these dirs are writable by the `users` group.
RUN echo "auth requisite pam_deny.so" >> /etc/pam.d/su \
    && sed -i.bak -e 's/^%admin/#%admin/' /etc/sudoers \
    && sed -i.bak -e 's/^%sudo/#%sudo/' /etc/sudoers \
    && useradd -m -s /bin/bash -N -u $NB_UID $NB_USER \
    && groupmod -g 1000 users \
    && mkdir -p $CONDA_DIR \
    && chown $NB_USER:$NB_GID $CONDA_DIR  \
    && chmod g+w /etc/passwd
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# INSTALL CONDA
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
  && /bin/bash Miniforge3-Linux-x86_64.sh -f -b -p "${CONDA_DIR:?}" \
  && rm Miniforge3-Linux-x86_64.sh \
  && $CONDA_DIR/bin/conda config --system --set auto_update_conda false \
  && $CONDA_DIR/bin/conda config --system --set show_channel_urls true \
  && $CONDA_DIR/bin/conda config --system --set auto_activate_base false \
  && $CONDA_DIR/bin/conda config --system --prepend channels conda-forge \
  && $CONDA_DIR/bin/conda config --system --add allowlist_channels conda-forge \
  && $CONDA_DIR/bin/conda config --system --set channel_priority strict \
  && $CONDA_DIR/bin/conda clean --all -f -y   \
  && find $CONDA_DIR -follow -type f -name '*.pyc' -delete \
  && /usr/local/bin/fix-permissions $HOME  \
  && /usr/local/bin/fix-permissions $CONDA_DIR
ENV PATH=$CONDA_DIR/bin:$PATH

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# INSTALL THE COCOA CONDA ENVIRONMENT (conda-lock yml downloaded from the
# COCOA_TAG git tag, so the environment matches the Cocoa version installed)

COPY init_cocoa.sh $HOME/init_cocoa.sh

RUN conda create -n lockenv -c conda-forge python=3.10 conda-lock=2.* wget --quiet --yes \
  && conda run --no-capture-output -n lockenv wget -O $HOME/cocoapy311-linux.yml \
     "https://raw.githubusercontent.com/CosmoLike/cocoa/refs/tags/${COCOA_TAG}/cocoapy311-linux.yml" \
  && conda run --no-capture-output -n lockenv conda-lock install --name cocoa $HOME/cocoapy311-linux.yml \
  && conda run --no-capture-output -n cocoa /bin/bash $HOME/init_cocoa.sh \
  && find $CONDA_DIR -follow -type f -name '*.pyc' -delete \
  && rm -rf $CONDA_DIR/share/man/* \
  && rm -rf $CONDA_DIR/envs/cocoa/man/* \
  && /usr/local/bin/fix-permissions $HOME \
  && /usr/local/bin/fix-permissions $CONDA_DIR

ENV PATH=$CONDA_DIR/envs/cocoa/bin:$PATH
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
RUN mkdir -p $HOME/.jupyter

RUN /usr/local/bin/fix-permissions $CONDA_DIR
RUN /usr/local/bin/fix-permissions $HOME
RUN /usr/local/bin/fix-permissions $HOME/.jupyter

RUN chown -R $NB_USER:$NB_GID $HOME
RUN chown -R $NB_USER:$NB_GID $CONDA_DIR
RUN chown -R $NB_USER:$NB_GID $HOME/.jupyter

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# VM: solves an important use: ls -l `which sh` -> bash not dash
RUN echo "dash dash/sh boolean false" | debconf-set-selections

COPY jupyter_notebook_config.py $HOME/.jupyter/jupyter_notebook_config.py
RUN chown -R $NB_USER:$NB_GID $HOME/.jupyter/jupyter_notebook_config.py
COPY bashrc $HOME/bash.bashrc

RUN chmod +x /usr/bin/tini \
  && chmod a+rwx $HOME/bash.bashrc \
  && mkdir -p $HOME/host/ && chmod -R a+rwx $HOME/host/ \
  && rm -rf ~/.cache/pip

EXPOSE 8888

ENTRYPOINT ["tini", "-g", "--"]
CMD ["/bin/bash", "-c", "source $HOME/bash.bashrc; /bin/bash"]

USER $NB_UID
WORKDIR $HOME
VOLUME ["$HOME/host/"]

# HOW TO BUILD AND PUBLISH THIS IMAGE
# docker build . --build-arg COCOA_TAG=v4.11.8 -t vivianmiranda/whovian-cocoa
# docker push vivianmiranda/whovian-cocoa
#
# HOW TO PUBLISH THE :thin VARIANT (image with Cocoa preinstalled/precompiled)
# run the container, then install and compile cocoa inside it
# docker ps -a                                          # check the container id
# docker commit <containerid> vivianmiranda/whovian-cocoa:thin
# docker push vivianmiranda/whovian-cocoa:thin
