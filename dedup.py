import sys
from sequniq import records


for record in records.uniqseqs(sys.stdin):
    print "%s\n%s\n+\n%s\n%s\n%s\n+\n%s" % record
