#!/usr/bin/env python
import subprocess

verbose=True
command = f'java -jar stilts/stilts.jar -h'
proc = subprocess.run(command, shell=True,
                        capture_output=(not verbose),
                        encoding="utf-8", check=True)