
  
import os
import re
import subprocess

from pysam import AlignmentFile #pylint: disable=no-name-in-module
from pysam import VariantFile #pylint: disable=no-name-in-module
import pandas as pd
import requests
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-vcf", "--vcf", type=str, required=True, help="input vcf")
parser.add_argument("-cns", "--cns", type=str, required=True, help="input cns")
parser.add_argument("-o", "--output", required=True, type=str, help="output directory")

def chromosome(chrom):
    chrs={"chr1":1,"chr2":2,"chr3":3,"chr4":4,"chr5":5,"chr6":6,"chr7":7,"chr8":8,"chr9":9,"chr10":10,"chr11":11,"chr12":12,"chr13":13,"chr14":14,"chr15":15,"chr16":16,"chr17":17,"chr18":18,"chr19":19,"chr20":20,"chr21":21,"chr22":22,"chrX":"X","chrY":"Y","NC_012920.1":"MT"}
    ch=chrs[chrom]
    return str(ch)

class SV:
    def __init__(self, rec):
        """ Create SV object for each pysam variant record."""
        self.chr1 = chromosome(rec.chrom)
        self.pos1 = int(rec.start) + 1
        self.pos2= int(rec.stop)
        self.type = rec.info["SVTYPE"]
        if self.type  == "BND":
            self.chr2=chromosome(rec.info["CHR2"])
        else:
            self.chr2=self.chr1

        self.name = "%s(%s:%s-%s:%s)" % (
            self.type,
            self.chr1,
            self.pos1,
            self.chr2,
            self.pos2,
        )

class Sample:
    """Creation of sample object."""
    def __init__(
        self,
        cns,
        merged_vcf,
        out_dir,
    ):
        """Creation of sample object."""
        # pylint: disable=line-too-long
        self.cns = cns
        self.out_dir = os.path.join(out_dir)
        self.vcf = str(merged_vcf)
        self.segs = os.path.join(self.out_dir, "segs.csv")
        self.circos_out = os.path.join(self.out_dir, "circos.png")
        self.svs = {}

        if os.path.isfile(self.vcf):
            self.svs = self.get_svs()

    def get_svs(self):
        """Get svs objects from sample vcf."""
        sv_list = [ SV(rec) for rec in VariantFile(self.vcf).fetch()]
        sv_name_dict = {sv.name: sv for sv in sv_list}
        return sv_name_dict

    def runcircos(self):
        """Set arguments and run script for creating circos plot."""
        pd.read_csv(self.cns, sep="\t")[
            ["chromosome", "start", "end", "tcn"]
        ].rename({"chromosome": "chrm", "tcn": "cns"}, axis=1).to_csv(
            self.segs, index=None
        )

        passed_svs = [
            sv
            for sv in self.svs.values()
        ]
        circos_sv_file = os.path.join(
            self.out_dir, "circos_svs.tsv"
        )
        circos_df = pd.DataFrame(
            [
                ("chr" + sv.chr1, sv.pos1, sv.pos1, "chr" + sv.chr2, sv.pos2, sv.pos2)
                for sv in passed_svs
            ],
            columns=[
                "Chromosome",
                "chromStart",
                "chromEnd",
                "Chromosome.1",
                "chromStart.1",
                "chromEnd.1",
            ],
        )
        circos_df.to_csv(circos_sv_file, index=None)


if __name__ == "__main__":
    args = parser.parse_args()
    sampl= Sample(args.cns,args.vcf,args.o)
    sampl.runcircos()