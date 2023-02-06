"""Tests for SPICE init."""

import unittest
import subprocess as sp
import socket
import time

def start_mongodb(port: int = 27017) -> str:
    """Start a local MongoDB instance."""
    cmd = f"mongod --dbpath data/testdb --port {port}"
    sp.run(cmd.split(), check=True)
    time.sleep(1)
    host = socket.gethostbyname(socket.gethostname())
    return host

