#!/bin/bash
# Convert single-particle REC.slcio to h5 with the same converter as nnqq training
# (LCIO2ak2_brems.py). Fixed-jet samples are evaluation-only; do not use pandoraFixed here.
#
# Usage:
#   ./makeak_single_particle.sh <particle>
#
# particle : gamma | neutron | kaon_long | pion | pi+ | pi- | muon | mu+ | mu- | electron | e+ | e-
#            kaon is accepted as an alias of kaon_long
#
# Input:
#   gpfs/pfa/simulation/generation/data/<particle>/reco/*_reco_REC.slcio
# Output:
#   gpfs/pfa/data/skimmed/pandora/singleParticle/<particle>/awkd/concat/<particle>_<i>.h5
#   (neutron_1000_reco_REC.slcio -> neutron_1.h5; neutron_0_reco_REC.slcio -> neutron_0.h5)
#
# Example:
#   ./makeak_single_particle.sh neutron

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "usage: $0 <particle>" >&2
  echo "  particle: gamma|neutron|kaon_long|pion|pi+|pi-|muon|mu+|mu-|electron|e+|e-" >&2
  exit 1
fi
particle=$1
if [ "$particle" = "kaon" ]; then
  particle=kaon_long
fi

scriptdir=$(cd "$(dirname "$0")" && pwd)
simroot=$(cd "${scriptdir}/.." && pwd)
cd "${simroot}"

converter=${simroot}/LCIO2ak2_brems.py
recodir=$(cd "${simroot}/../gpfs/pfa/simulation/generation/data/${particle}/reco" && pwd)
datadir=$(cd "${simroot}/../gpfs/pfa/data/skimmed/pandora" && pwd)/singleParticle/${particle}
jobdir=${simroot}/job/skimmed/pandora/singleParticle/${particle}/concat

if [ ! -f "${converter}" ]; then
  echo "converter not found: ${converter}" >&2
  exit 1
fi
if [ ! -d "${recodir}" ]; then
  echo "reco dir not found: ${recodir}" >&2
  exit 1
fi

mkdir -p "${datadir}/log/concat" "${datadir}/awkd/concat" "${jobdir}"

shopt -s nullglob
files=("${recodir}"/*_reco_REC.slcio)
if [ ${#files[@]} -eq 0 ]; then
  echo "no *_reco_REC.slcio in ${recodir}" >&2
  exit 1
fi

nsubmit=0
nskip=0
for file in "${files[@]}"; do
  filename=$(basename "${file}" .slcio)
  name=${filename%000_reco_REC}
  if [ "${name}" = "${filename}" ]; then
    name=${filename%_reco_REC}
  fi
  outh5=${datadir}/awkd/concat/${name}.h5
  logfile=${datadir}/log/concat/${name}.log

  echo "${file}"
  if [ -e "${outh5}" ]; then
    echo "   ${outh5} already exist"
    nskip=$((nskip + 1))
    continue
  fi

  event_count=$(lcio_event_counter "${file}" 2>/dev/null | tr -dc '0-9' || true)
  if [ -z "${event_count}" ] || [ "${event_count}" -eq 0 ]; then
    echo "   skip (no events): ${file}" >&2
    continue
  fi

  echo "   ${name}.h5  (${event_count} events)"
  bsub -q s -o "${jobdir}/output.%J" -e "${jobdir}/errors.%J" \
    "python ${converter} ${file} ${outh5} ${event_count} 0 > ${logfile}"
  nsubmit=$((nsubmit + 1))
done

echo "submitted ${nsubmit} jobs, skipped ${nskip} existing h5 (${#files[@]} REC files)"
