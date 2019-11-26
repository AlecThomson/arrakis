#!/usr/bin/env python
import subprocess
from pathlib import Path
import os

verbose=True
stiltspath = Path(os.path.realpath(__file__)).parent.parent/"stilts"/"stilts.jar"
command = ['java','-jar', stiltspath,'-h']
proc = subprocess.run(command,
                        capture_output=(not verbose),
                        encoding="utf-8", check=True)