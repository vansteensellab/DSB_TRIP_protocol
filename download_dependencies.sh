#!/usr/bin/env bash

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
mkdir -p $SCRIPT_DIR/lib
git clone git@github.com:maxwshen/inDelphi-model.git $SCRIPT_DIR/lib/inDelphi-model
git clone git@github.com:felicityallen/SelfTarget.git $SCRIPT_DIR/lib/SelfTarget
