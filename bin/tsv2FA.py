#!/usr/bin/env python
import sys
import os
def main():
    In = sys.argv[1]
    Out = sys.argv[2]
    f = open(In,"rb")
    data = f.readlines()
    f.close()
    f = open(Out,"w")
    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        f.writelines([">" + tmp[1] + os.linesep])
        f.writelines([tmp[2] + os.linesep])
    f.close()
if __name__ == "__main__":
    main()
