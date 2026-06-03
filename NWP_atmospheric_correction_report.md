# NWP Atmospheric Correction Support Report (Repository Evidence)

## Scope
This report summarizes repository evidence for which NWP products are implemented and used for atmospheric correction in ECOSTRESS Collection 3 L2 LSTE, with focus on MERRA-2, GEOS/GEOS-5/GEOS-5 FP, NCEP, and ECMWF.

## 1) Implemented NWP products and usage in atmospheric correction

### Implemented in code
- **MERRA-2** is implemented and used in an operational branch:
  - Selection/dispatch: `tes_read_interp_atmos_merra(...)` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1387-1401`)
  - MERRA2 file matching: `*MERRA2*inst6*Np*...nc4` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1317`, `1391`)
  - MERRA2 validation in filename: (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1354`, `1428`)

- **GEOS (GEOS5-style / GEOS FP filename patterns)** is implemented and used:
  - Selection/dispatch: `tes_read_interp_atmos_geos(...)` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1435-1448`)
  - GEOS file matching pattern: `GEOS*3d_asm_Np*...nc4` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1064`)
  - GEOS source metadata value written as `"GEOS5"` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1476-1478`)

- **NCEP** is implemented and used:
  - Selection/dispatch: `tes_read_interp_atmos_ncep(...)` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1486-1501`)
  - NCEP filename pattern: `gdas1*...` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1858`)
  - NCEP GRIB parsing pipeline in `read_atmos_ncep(...)` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1986-2216`)
  - Build requires ecCodes/GRIB API (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/Makefile:181-206`)

- **ECMWF** is listed but not implemented operationally:
  - Explicit error: `ECMWF processing is not implemented yet.` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1684`)
  - Placeholder TODO for ECMWF read function (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1693`)

### Atmospheric correction linkage
NWP ingest is performed directly before RTTOV profile preparation/execution, so selected NWP feeds atmospheric correction: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1374-1730`.

## 2) Exact code paths for NWP source detection and selection

1. Read `NWP_DIR` from run config `StaticAncillaryFileGroup`:
   - `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:574-583`

2. Read source keys from run parameters group `NWP` (`ecmwf`, `geos`, `merra`, `ncep`):
   - `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:855-891`
   - `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/OSP/PgeRunParameters.xml:79-86`

3. Detect source by substring search in `NWP_DIR`:
   - `strstr(NWP_dir, nwp_key_...)` at `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1379-1382`

4. Dispatch branches:
   - MERRA: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1387-1434`
   - GEOS: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1435-1485`
   - NCEP: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1486-1678`
   - ECMWF fail path: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1679-1720`
   - Unknown source hard fail: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1723-1725`

## 3) Evidence by product: documented support vs implemented code vs operational/default configuration

### MERRA / MERRA-2
- **Documented support**:
  - NWP selectable keys include MERRA (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/README.md:197`)
  - Spec describes MERRA branch (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/ECOSTRESS Collection 3 TES Algorithm Specification.md:218-223`)
- **Implemented code**: yes (MERRA2 file naming/ingest paths in `tes_util.c` and branch in `tes_main.c`; citations above).
- **Sample/default config**: not the sample default (sample points to GEOS5; see Section 4).

### GEOS / GEOS5 / GEOS-5 FP
- **Documented support**:
  - User guide says atmospheric correction uses GEOS5 profiles (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/documentation/ECOL2_User_Guide_V3.md:242`)
  - Input ancillary table lists **GEOS5-FP** (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/documentation/ECOL2_User_Guide_V3.md:261`)
  - Cloud algorithm text references GEOS5-FP (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/documentation/ECOL2_User_Guide_V3.md:418`)
- **Implemented code**:
  - GEOS branch and GEOS file pattern matching (`src/tes_main.c`, `src/tes_util.c`; above)
  - Default ancillary pointer string uses GEOS FP style filename (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:317`)
- **GEOS-5 FP data products and layers actually read**:
  - Product files searched/read: GEOS 3-hourly `3d_asm_Np` NetCDF files matching `GEOS*3d_asm_Np*YYYYMMDD_hh00*.nc4` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1063-1065`).
  - Two time-adjacent files are read per scene (before and after granule time), then interpolated (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1139-1176`, `1202-1212`).
  - Data layers read from each GEOS file: `lev`, `lat`, `lon`, `T`, `QV`, `PS` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:1023-1035`).
  - `tcw` is not directly read from GEOS in this path; if missing, TCW is derived later from humidity/pressure (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1739-1758`).
- **Sample/default config**:
  - Sample run config uses GEOS5 path (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/config/sample_RunConfig.xml:11`)

### NCEP
- **Documented support**:
  - README/spec list NCEP as selectable NWP source (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/README.md:197`; `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/ECOSTRESS Collection 3 TES Algorithm Specification.md:218-223`)
- **Implemented code**: yes (NCEP branch and GRIB ingest/conversion, citations above).
- **Sample/default config**: not sample default.

### ECMWF
- **Documented support**:
  - Listed as selectable key in runtime parameters and docs (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/OSP/PgeRunParameters.xml:82`; `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/README.md:197`)
- **Implemented code**: no (explicit not implemented error and TODO in `tes_main.c`).
- **Operational/default config**: no.

## 4) Sample/default run configuration
- Sample run config sets:
  - `NWP_DIR = /NWP/GEOS5` (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/config/sample_RunConfig.xml:11`)
- Runtime NWP keys are:
  - `ecmwf=ECMWF`, `geos=GEOS`, `merra=MERRA`, `ncep=NCEP`
  - (`/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/OSP/PgeRunParameters.xml:82-85`)

Therefore, with provided sample config, source detection resolves to the GEOS branch.

## 5) Does this repository prove actual ECOSTRESS production uses GEOS-5 FP?

### What this repository does show
- Code supports multiple NWP sources and selects by runtime-configured `NWP_DIR` + key matching.
- Sample/default example config points to GEOS5.
- User documentation in-repo describes GEOS5/GEOS5-FP atmospheric inputs.

### What this repository does not prove by itself
- It does **not** by itself prove live operational production always uses GEOS-5 FP, because actual production selection is runtime-configurable and external deployment/orchestration settings are outside this codebase.
  - Runtime-configurable selection evidence: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:574-583`, `1379-1382`.

## 6) Metadata handling evidence relevant to NWP source/file traceability
- `NWPSource` is set in branch-specific logic:
  - MERRA2: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1425-1427`
  - GEOS5: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1476-1478`
  - NCEP: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1663-1665`
  - ECMWF: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1711-1713`
- `AncillaryInputPointer` is written to output metadata and updated with GEOS filenames in GEOS branch:
  - Initialization default GEOS FP-style string: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:317`
  - GEOS branch update: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:1454-1455`
  - Metadata writeout: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:4341-4347`
- `AncillaryNWP` metadata is populated from internal tracked NWP file name(s):
  - `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_main.c:4542-4543`
  - `get_nwp_name()` backing state: `/tmp/workspace/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/src/tes_util.c:54-57`

## 7) Conclusions
1. **Implemented NWP ingestion for atmospheric correction:** MERRA-2, GEOS, and NCEP are implemented and wired into RTTOV preprocessing.  
2. **ECMWF:** listed as a selectable source but currently not operationally implemented (explicit error path).  
3. **Sample/default config in repository:** GEOS5.  
4. **GEOS-5 FP for production claim:** repository evidence strongly supports GEOS5/GEOS5-FP intent and sample configuration, but does not alone prove all real production runs use GEOS-5 FP; that requires external operational evidence.

## Confidence and limitations
- **Confidence: High** for code-level implementation status and runtime selection behavior.
- **Confidence: High** that sample/default config is GEOS5.
- **Limitation:** no external production runtime/deployment configuration or operations records are included in this repository.
