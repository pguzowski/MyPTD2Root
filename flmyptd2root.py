#!/usr/bin/env python
import os
import argparse
import subprocess
from contextlib import contextmanager
import shutil
import tempfile

@contextmanager
def tmpdir():
    temp_dir = tempfile.mkdtemp()
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir)

def getConversionScript():
    """ derive path to conversion script from this
    """
    thisScript = os.path.realpath(__file__)
    return os.path.join(os.path.dirname(thisScript), "myp2r.conf")


def flbrio2root():
    """implementation of Brio->Root conversion
    """
    parser = argparse.ArgumentParser(description="Convert flreconstruct BRIO format to ROOT")
    parser.add_argument("-i", type=str, metavar="<infile>", required=True, help="input BRIO file")
    parser.add_argument("-o", type=str, metavar="<outfile>", required=True, help="output ROOT file")
    args = parser.parse_args()

    fullinpath = os.path.realpath(args.i)
    fulloutpath = os.path.realpath(args.o)

    with tmpdir() as temp_dir:

        try:
            subprocess.call(["flreconstruct", "-i", fullinpath, "-p", getConversionScript()],cwd=temp_dir)
            shutil.move(temp_dir+'/testptd.root', fulloutpath)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                print("error: could not locate flreconstruct program")
            else:
                # Something else went wrong while trying to run `flreconstruct`
                raise


if __name__ == '__main__':
    flbrio2root()
