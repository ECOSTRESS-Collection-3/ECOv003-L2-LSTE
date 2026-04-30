# Root Makefile — delegates all targets to src/Makefile
# Executable is built to bin/L2_PGE
#
# Targets:
#   all     build the L2_PGE executable (default)
#   clean   remove compiled objects and executable
#   environment create mamba env if missing
#   install copy executable to PREFIX/bin (default PREFIX=/usr/local)

PREFIX ?= /usr/local
ENV_NAME ?= ECOv003-L2-LSTE
MAMBA ?= mamba
ENV_PACKAGES := -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config

all:
	$(MAKE) -C src DEFAULT_ENV_NAME="$(ENV_NAME)" all

clean:
	$(MAKE) -C src clean

environment:
	@$(MAMBA) --version >/dev/null 2>&1 || { \
		echo "Error: '$(MAMBA)' not found in PATH."; \
		echo "Install mamba (or micromamba) and retry."; \
		exit 1; \
	}
	@if $(MAMBA) env list 2>/dev/null | awk '{print $$1}' | grep -Fxq "$(ENV_NAME)"; then \
		echo "Environment '$(ENV_NAME)' already exists; skipping create."; \
	else \
		echo "Creating environment '$(ENV_NAME)'..."; \
		$(MAMBA) create -n "$(ENV_NAME)" $(ENV_PACKAGES) -y; \
	fi

install: environment all
	$(MAKE) -C src DEFAULT_ENV_NAME="$(ENV_NAME)" install PREFIX="$(PREFIX)" DESTDIR="$(DESTDIR)"

.PHONY: all clean environment install
