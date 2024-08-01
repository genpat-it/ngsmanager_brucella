#!/bin/bash
set -euo pipefail

LOCAL_REGISTRY=nexus-prod.izs.intra:9091

grep -r "container " ../[sm]*/*.nf | sed 's/^.*container //g' | sed "s/['\"]//g" | sed "s/..LOCAL_REGISTRY./${LOCAL_REGISTRY}/g" | sort -u