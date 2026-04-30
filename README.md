# ECOSTRESS Level 2 Surface Temperature

This will be the repository for the ECOsystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) collection 3 level 2 surface temperature data product algorithm.

The ECOSTRESS collection 3 level 2 surface temperature data product is the pre-cursor to the [Surface Biology and Geology (SBG) collection 1 level 2 surface temperature data product algorithm](https://github.com/sbg-tir/SBG-TIR-L2-LSTE).

## Repository Structure

```
ECOv003-L2-LSTE/
├── Makefile          # top-level build entry point
├── README.md
├── LICENSE
├── bootstrap.sh      # auto-installs dependencies and builds
├── include/          # C header files
├── src/              # C and Fortran source files
├── scripts/          # run_ecostress_fwd_*.sh helper scripts
├── config/           # sample configuration XML files
├── OSP/              # look-up tables and coefficient files
├── build/            # compiled object files (generated, not committed)
└── bin/              # compiled executable (generated, not committed)
```

## Building

### ✓ Recommended: One-Command Install (using Mamba)

This is the **easiest and most portable approach** that works on macOS, Linux, and Windows:

```bash
# Create the mamba environment if missing
make environment

# Activate the environment
mamba activate ECOv003-L2-LSTE

# Build and install from the repo root
make install
```

The executable is built to `bin/L2_PGE` and installed to `/usr/local/bin` by default.

**Quick one-liner:**
```bash
make environment && mamba activate ECOv003-L2-LSTE && make install
```

For a user-local install (no sudo):

```bash
make environment && mamba activate ECOv003-L2-LSTE && make install PREFIX=$HOME/.local
```

### Manual Build

If you already have all dependencies installed:

```bash
make
```

### Build Commands Reference

```bash
# Create/update dependencies environment
make environment

# Clean build artifacts
make clean

# Regular build only (no install)
make

# Build and install (recommended)
make install

# Debug build (with symbols, no optimization)
DEBUG=1 make
```

### Troubleshooting

**"HDF4 not found"**
- **Solution:** Use mamba—HDF4 is not available in Homebrew
- Create environment: `make environment`

**"libxml2 not found"** or **"expat not found"**
- **Solution:** Use mamba (see above)—these packages don't always provide pkg-config files
- Mamba handles this automatically

**"pkg-config: command not found"**
- macOS: `brew install pkg-config` or use mamba (included in environment)
- Ubuntu: `sudo apt-get install pkg-config`
- CentOS: `sudo yum install pkgconfig`

**Build fails with "format specifies type" warnings**
- These are harmless; the executable is still created and works
- To suppress: add `-Wno-format` to `ADD_CFLAGS` in `src/Makefile`

**On macOS with M1/M2 (ARM):**
- Mamba automatically installs ARM-native binaries—no special action needed

**Still having issues?**
- Clean and rebuild: `make clean && make`
- Verify environment: `mamba list -n ECOv003-L2-LSTE | grep -E "hdf|eccodes|libxml"`
- Check pkg-config: `pkg-config --list-all | grep hdf`
