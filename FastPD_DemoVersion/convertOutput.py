import argparse
import uu
import binascii

def main():
    parser= argparse.ArgumentParser(description="Get file")
    parser.add_argument("output")
    args= parser.parse_args()

#    uu.decode(args.output, args.output+".txt")

    outfile= open(args.output+".txt", "w")
    f= open(args.output, "rb")
    try:
        byte=f.read(1)
        while byte != "":
#            num= 
            num= binascii.b2a_uu(byte)
            outfile.write(num+"\n")
            byte= f.read(1)
    finally:
        f.close()  
 
#    with open(args.output, "rb") as output:
#        data= output.read(1)
#        text= data.decode('utf-8')
#    for byte in text:
#         print byte
#         outfile.write(byte+"\n")

    outfile.close()
#    output.close()

main()
