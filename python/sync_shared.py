"""Copy the files shared with the R package (inst/shared: labels, French help, reference values) into the
Python package before building or testing. Usage: python3 sync_shared.py"""
import os, shutil
here = os.path.dirname(os.path.abspath(__file__))
src, dst = os.path.join(here, "..", "inst", "shared"), os.path.join(here, "canofac", "shared")
shutil.rmtree(dst, ignore_errors=True)
shutil.copytree(src, dst)
print("copied", src, "->", dst)
