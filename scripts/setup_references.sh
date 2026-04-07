#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${FINALEME_DATA_DIR:-${ROOT_DIR}/data}"
ZENODO_URL="https://zenodo.org/records/19392525/files"
JAR_PATH="${ROOT_DIR}/target/FinaleMe-0.61-jar-with-dependencies.jar"

OS=""
ARCH=""
UCSC_PLATFORM=""

log() {
  printf '[setup] %s\n' "$*"
}

warn() {
  printf '[setup][warn] %s\n' "$*" >&2
}

die() {
  printf '[setup][error] %s\n' "$*" >&2
  exit 1
}

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

is_maven_3_8_or_newer() {
  local version="$1"
  local major minor
  major="$(printf '%s' "${version}" | awk -F. '{print $1}')"
  minor="$(printf '%s' "${version}" | awk -F. '{print $2}')"
  if [[ -z "${major}" || -z "${minor}" ]]; then
    return 1
  fi
  if (( major > 3 )); then
    return 0
  fi
  if (( major == 3 && minor >= 8 )); then
    return 0
  fi
  return 1
}

detect_platform() {
  OS="$(uname -s)"
  ARCH="$(uname -m)"

  case "${OS}" in
    Darwin)
      case "${ARCH}" in
        arm64) UCSC_PLATFORM="macOSX.arm64" ;;
        x86_64) UCSC_PLATFORM="macOSX.x86_64" ;;
        *) die "Unsupported macOS architecture: ${ARCH}" ;;
      esac
      ;;
    Linux)
      case "${ARCH}" in
        x86_64) UCSC_PLATFORM="linux.x86_64" ;;
        aarch64)
          UCSC_PLATFORM="linux.x86_64"
          warn "Linux arm64 detected; UCSC may not provide native arm64 binaries. Using linux.x86_64 URLs."
          ;;
        *) die "Unsupported Linux architecture: ${ARCH}" ;;
      esac
      ;;
    *)
      die "Unsupported OS: ${OS}" ;;
  esac

  log "Detected platform: OS=${OS}, ARCH=${ARCH}, UCSC_PLATFORM=${UCSC_PLATFORM}"
}

check_java() {
  if ! command_exists java; then
    die "Missing required dependency: java. Install Oracle JDK 21: https://www.oracle.com/java/technologies/downloads/#java21"
  fi

  local java_version major
  java_version="$(java -version 2>&1 | awk -F'"' '/version/ {print $2; exit}')"
  if [[ -z "${java_version}" ]]; then
    java_version="$(java -version 2>&1 | head -n 1)"
  fi

  major="$(printf '%s' "${java_version}" | awk -F. '{if ($1=="1") print $2; else print $1}' | sed 's/[^0-9].*$//')"
  if [[ -z "${major}" ]]; then
    die "Unable to parse Java version from: $(java -version 2>&1 | head -n 1)"
  fi
  if (( major < 21 )); then
    die "Java ${major} detected; Java 21+ required. Install Oracle JDK 21: https://www.oracle.com/java/technologies/downloads/#java21"
  fi

  log "Java version OK: ${java_version}"
}

check_deps() {
  log "Checking dependencies..."
  check_java

  if ! command_exists git; then
    die "Missing required dependency: git"
  fi

  if ! command_exists curl && ! command_exists wget; then
    die "Missing download client: install curl or wget"
  fi

  if command_exists mvn; then
    local mvn_version
    mvn_version="$(mvn -v 2>/dev/null | awk '/Apache Maven/ {print $3; exit}')"
    log "Maven detected: $(mvn -v 2>/dev/null | head -n 1)"
    if [[ -n "${mvn_version}" ]] && ! is_maven_3_8_or_newer "${mvn_version}"; then
      warn "Maven ${mvn_version} detected; Maven 3.8+ is recommended for source builds."
    fi
  else
    warn "Maven not found. Needed only for 'build' or full setup from source. Install Maven 3.8+ or use a prebuilt JAR from Releases."
  fi

  if ! command_exists samtools; then
    warn "samtools not found (needed for Step 1 FASTA/chrom-size helper workflows)."
  fi
  if ! command_exists bedtools; then
    warn "bedtools not found (needed for Step 4 helper workflows)."
  fi
  if ! command_exists bgzip; then
    warn "bgzip not found (needed for fragment bgzip workflows)."
  fi
  if ! command_exists tabix; then
    warn "tabix not found (needed for fragment tabix workflows)."
  fi

  if ! command_exists bedGraphToBigWig; then
    warn "bedGraphToBigWig not found (needed for FinaleMe -bwOutput in Step 3)."
    warn "Download hint: https://hgdownload.soe.ucsc.edu/admin/exe/${UCSC_PLATFORM}/bedGraphToBigWig"
  fi

  if ! command_exists twoBitInfo; then
    warn "twoBitInfo not found; chrom.sizes will be downloaded from UCSC instead of generated from .2bit."
    warn "Download hint: https://hgdownload.soe.ucsc.edu/admin/exe/${UCSC_PLATFORM}/twoBitInfo"
  fi

  warn "Java install recommendation: Oracle JDK 21 from https://www.oracle.com/java/technologies/downloads/#java21"
  warn "Pick your platform (macOS arm64/x64, Linux x64/aarch64). On clusters without root, use the tarball and set JAVA_HOME/PATH."
}

download_file() {
  local url="$1"
  local dest="$2"
  local tmp="${dest}.tmp"

  if [[ -s "${dest}" ]]; then
    log "Skipping existing: ${dest}"
    return
  fi

  mkdir -p "$(dirname "${dest}")"
  rm -f "${tmp}"

  if command_exists curl; then
    log "Downloading with curl: ${url}"
    curl -fSL -o "${tmp}" "${url}"
  elif command_exists wget; then
    log "Downloading with wget: ${url}"
    wget -q -O "${tmp}" "${url}"
  else
    die "Neither curl nor wget is available"
  fi

  mv "${tmp}" "${dest}"
  local bytes
  bytes="$(wc -c < "${dest}" | tr -d '[:space:]')"
  log "Saved ${dest} (${bytes} bytes)"
}

build_finaleme() {
  if [[ -s "${JAR_PATH}" ]]; then
    log "Skipping build; JAR already exists: ${JAR_PATH}"
    return
  fi

  if ! command_exists mvn; then
    die "Maven is required to build FinaleMe. Install Maven 3.8+ or download a prebuilt JAR from GitHub Releases."
  fi
  local mvn_version
  mvn_version="$(mvn -v 2>/dev/null | awk '/Apache Maven/ {print $3; exit}')"
  if [[ -z "${mvn_version}" ]] || ! is_maven_3_8_or_newer "${mvn_version}"; then
    die "Maven 3.8+ is required for building FinaleMe. Detected: ${mvn_version:-unknown}"
  fi

  log "Syncing vendored Maven repository..."
  "${ROOT_DIR}/scripts/sync-vendored-repo.sh"

  log "Building FinaleMe JAR..."
  mvn -f "${ROOT_DIR}/pom.xml" clean package -DskipTests

  if [[ ! -s "${JAR_PATH}" ]]; then
    die "Build finished but JAR not found: ${JAR_PATH}"
  fi
  log "Build complete: ${JAR_PATH}"
}

setup_genomes() {
  log "Setting up genome .2bit files in ${DATA_DIR}"
  download_file "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit" "${DATA_DIR}/hg19.2bit"
  download_file "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit" "${DATA_DIR}/hg38.2bit"
}

setup_chrom_sizes() {
  log "Setting up chromosome size files"

  local hg19_sizes="${DATA_DIR}/hg19.chrom.sizes"
  local hg38_sizes="${DATA_DIR}/hg38.chrom.sizes"

  if command_exists twoBitInfo; then
    setup_genomes

    if [[ -s "${hg19_sizes}" ]]; then
      log "Skipping existing: ${hg19_sizes}"
    else
      twoBitInfo "${DATA_DIR}/hg19.2bit" "${hg19_sizes}"
      log "Generated ${hg19_sizes} via twoBitInfo"
    fi

    if [[ -s "${hg38_sizes}" ]]; then
      log "Skipping existing: ${hg38_sizes}"
    else
      twoBitInfo "${DATA_DIR}/hg38.2bit" "${hg38_sizes}"
      log "Generated ${hg38_sizes} via twoBitInfo"
    fi
  else
    warn "twoBitInfo not available; downloading prebuilt chrom.sizes files"
    download_file "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes" "${hg19_sizes}"
    download_file "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" "${hg38_sizes}"
  fi
}

setup_cpg_files() {
  log "Downloading CpG motif and index files from Zenodo"
  download_file "${ZENODO_URL}/CG_motif.hg19.common_chr.pos_only.bedgraph.gz?download=1" "${DATA_DIR}/CG_motif.hg19.common_chr.pos_only.bedgraph.gz"
  download_file "${ZENODO_URL}/CG_motif_seqkit.pos_only.hg38.bedgraph.gz?download=1" "${DATA_DIR}/CG_motif_seqkit.pos_only.hg38.bedgraph.gz"
  download_file "${ZENODO_URL}/CpG_index.hg19.bed.gz?download=1" "${DATA_DIR}/CpG_index.hg19.bed.gz"
  download_file "${ZENODO_URL}/CpG_index.hg19.bed.gz.csi?download=1" "${DATA_DIR}/CpG_index.hg19.bed.gz.csi"
  download_file "${ZENODO_URL}/CpG_index.hg38.bed.gz?download=1" "${DATA_DIR}/CpG_index.hg38.bed.gz"
  download_file "${ZENODO_URL}/CpG_index.hg38.bed.gz.csi?download=1" "${DATA_DIR}/CpG_index.hg38.bed.gz.csi"
}

setup_dark_regions() {
  log "Downloading dark-region BED files"
  download_file "${ZENODO_URL}/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed?download=1" \
    "${DATA_DIR}/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed"

  local hg38_blacklist_gz="${DATA_DIR}/hg38-blacklist.v2.bed.gz"
  local hg38_blacklist_bed="${DATA_DIR}/hg38-blacklist.v2.bed"
  download_file "${ZENODO_URL}/hg38-blacklist.v2.bed.gz?download=1" "${hg38_blacklist_gz}"

  if [[ -s "${hg38_blacklist_bed}" ]]; then
    log "Skipping existing: ${hg38_blacklist_bed}"
  else
    gzip -dc "${hg38_blacklist_gz}" > "${hg38_blacklist_bed}.tmp"
    mv "${hg38_blacklist_bed}.tmp" "${hg38_blacklist_bed}"
    log "Extracted ${hg38_blacklist_bed}"
  fi
}

setup_methylation_prior() {
  log "Downloading methylation prior bigWig files"
  download_file "${ZENODO_URL}/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw?download=1" "${DATA_DIR}/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw"
  download_file "${ZENODO_URL}/wgbs_buffyCoat_jensen2015GB.methy.hg38.bw?download=1" "${DATA_DIR}/wgbs_buffyCoat_jensen2015GB.methy.hg38.bw"
}

print_summary() {
  log "Reference data directory: ${DATA_DIR}"

  local files=(
    "hg19.2bit|reference genome (hg19)"
    "hg38.2bit|reference genome (hg38)"
    "hg19.chrom.sizes|chromosome sizes (hg19)"
    "hg38.chrom.sizes|chromosome sizes (hg38)"
    "CG_motif.hg19.common_chr.pos_only.bedgraph.gz|CpG motif coordinates for Step 1 (hg19)"
    "CG_motif_seqkit.pos_only.hg38.bedgraph.gz|CpG motif coordinates for Step 1 (hg38)"
    "CpG_index.hg19.bed.gz|CpG index for -cpgIndexFile (hg19)"
    "CpG_index.hg19.bed.gz.csi|tabix/CSI index for hg19 CpG index"
    "CpG_index.hg38.bed.gz|CpG index for -cpgIndexFile (hg38)"
    "CpG_index.hg38.bed.gz.csi|tabix/CSI index for hg38 CpG index"
    "wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed|dark regions for hg19"
    "hg38-blacklist.v2.bed|dark regions for hg38"
    "wgbs_buffyCoat_jensen2015GB.methy.hg19.bw|methylation prior (hg19)"
    "wgbs_buffyCoat_jensen2015GB.methy.hg38.bw|methylation prior (hg38)"
  )

  local missing=0
  printf '\n%-10s %-75s %s\n' "Status" "File" "Purpose"
  printf '%s\n' "---------------------------------------------------------------------------------------------------------------"
  for entry in "${files[@]}"; do
    local rel purpose path status
    rel="${entry%%|*}"
    purpose="${entry#*|}"
    path="${DATA_DIR}/${rel}"
    if [[ -s "${path}" ]]; then
      status="[ok]"
    else
      status="[missing]"
      missing=$((missing + 1))
    fi
    printf '%-10s %-75s %s\n' "${status}" "${path}" "${purpose}"
  done

  printf '\n'
  if (( missing > 0 )); then
    warn "${missing} file(s) are still missing."
  else
    log "All reference files are present."
  fi

  log "Note: wgbstools/UXM_deconv are only needed for optional custom atlas generation in Step 5."
}

usage() {
  cat <<USAGE
Usage: $(basename "$0") [command]

Commands:
  all         Run full setup (default)
  deps        Check system dependencies
  build       Build FinaleMe JAR
  genomes     Download hg19/hg38 .2bit files
  chromsizes  Generate/download chrom.sizes files
  cpg         Download CpG motif + index files
  darkregions Download dark-region BED files
  methylation Download methylation prior bigWig files
  summary     Print reference file summary
USAGE
}

run_all() {
  detect_platform
  check_deps
  build_finaleme
  setup_genomes
  setup_chrom_sizes
  setup_cpg_files
  setup_dark_regions
  setup_methylation_prior
  print_summary
}

main() {
  local cmd="${1:-all}"

  case "${cmd}" in
    all)
      run_all
      ;;
    deps)
      detect_platform
      check_deps
      ;;
    build)
      detect_platform
      check_deps
      build_finaleme
      ;;
    genomes)
      detect_platform
      setup_genomes
      ;;
    chromsizes)
      detect_platform
      setup_chrom_sizes
      ;;
    cpg)
      detect_platform
      setup_cpg_files
      ;;
    darkregions)
      detect_platform
      setup_dark_regions
      ;;
    methylation)
      detect_platform
      setup_methylation_prior
      ;;
    summary)
      detect_platform
      print_summary
      ;;
    -h|--help|help)
      usage
      ;;
    *)
      die "Unknown command: ${cmd}. Use --help for usage."
      ;;
  esac
}

main "$@"
