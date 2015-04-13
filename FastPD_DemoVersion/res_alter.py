from array import array
import struct
import numpy
import sys

out_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/results.bin', 'wb+')
out_file.write(struct.pack('=l', 1))
out_file.close()
print("ey yo we done")

