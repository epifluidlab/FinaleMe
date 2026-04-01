#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIB_DIR="${ROOT_DIR}/lib"
REPO_DIR="${ROOT_DIR}/lib-repo"

publish_artifact() {
  local group_id="$1"
  local artifact_id="$2"
  local version="$3"
  local source_jar="$4"

  local group_path
  group_path="$(printf '%s' "${group_id}" | tr '.' '/')"
  local target_dir="${REPO_DIR}/${group_path}/${artifact_id}/${version}"
  local target_jar="${target_dir}/${artifact_id}-${version}.jar"
  local target_pom="${target_dir}/${artifact_id}-${version}.pom"

  if [[ ! -f "${LIB_DIR}/${source_jar}" ]]; then
    echo "Missing source jar: ${LIB_DIR}/${source_jar}" >&2
    exit 1
  fi

  mkdir -p "${target_dir}"
  cp "${LIB_DIR}/${source_jar}" "${target_jar}"

  cat > "${target_pom}" <<EOF
<project xmlns="http://maven.apache.org/POM/4.0.0"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://maven.apache.org/POM/4.0.0 https://maven.apache.org/xsd/maven-4.0.0.xsd">
  <modelVersion>4.0.0</modelVersion>
  <groupId>${group_id}</groupId>
  <artifactId>${artifact_id}</artifactId>
  <version>${version}</version>
  <packaging>jar</packaging>
  <name>${artifact_id}</name>
</project>
EOF
}

publish_artifact "be.ac.ulg.montefiore.run.jahmm" "jahmm" "0.6.2" "jahmm-0.6.2.jar"
publish_artifact "edu.unc.genomics" "java-genomics-io" "1.0" "java-genomics-io.jar"
publish_artifact "org.broad.igv" "igv" "1.0" "igv.jar"
publish_artifact "ch.systemsx.cisd.hdf5" "sis-jhdf5-batteries_included" "1.0" "sis-jhdf5-batteries_included.jar"
publish_artifact "jdistlib" "jdistlib" "0.4.5" "jdistlib-0.4.5-bin.jar"

echo "Vendored Maven repository updated at: ${REPO_DIR}"
