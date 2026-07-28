#!/usr/bin/env bash
# Mamun raw deposit downloader (Materials Cloud xxdhq-f1v83, CC-BY 4.0).
set -e
cd "$(dirname "$0")/raw"
for f in references.tar.gz O.tar.gz C.tar.gz H.tar.gz N.tar.gz S.tar.gz; do
    echo "=== $f $(date +%H:%M:%S)"
    curl -sL --retry 3 \
        "https://archive.materialscloud.org/records/xxdhq-f1v83/files/$f?download=1" \
        -o "$f"
    du -h "$f"
done
echo ALL-DOWNLOADS-DONE
