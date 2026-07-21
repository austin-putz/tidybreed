# Detect available system memory in bytes (cross-OS, best-effort)

Returns available (not total) physical RAM in bytes on Linux, macOS
(Intel + Apple Silicon), and Windows, or `NA_real_` when detection fails
or the OS is unsupported. **Never errors** — callers fall back to a
conservative fixed budget when this returns `NA`.

## Usage

``` r
detect_available_memory()
```

## Value

Numeric scalar bytes, or `NA_real_`.

## Details

- Linux: `/proc/meminfo` `MemAvailable` (falls back to `MemFree`), in
  kB.

- macOS: `vm_stat` free + inactive + speculative pages x `hw.pagesize`.

- Windows: `wmic OS get FreePhysicalMemory` (kB), then a PowerShell
  fallback.
