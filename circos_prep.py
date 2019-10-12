
  
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
parser.add_argument("-ref", "--ref", type=str, required=True, help="reference bedfile")
parser.add_argument("-id", "--id", required=True, type=str, help="input id")
parser.add_argument("-o", "--output", required=True, type=str, help="output directory")

class SV:
    def __init__(self, rec):
        """ Create SV object for each pysam variant record."""
        self.chr1 = rec.contig
        self.pos1 = int(rec.start) + 1
        self.chr2, self.pos2 = re.sub(r"\[|\]|N", "", rec.alts[0]).split(
            ":"
        )  # can use chr2 and end (.stop)
        self.pos2 = int(self.pos2)
        info_keys = [k for k in rec.info.keys()]
        self.info_keys = info_keys
        if "SVTYPE" in info_keys:
            self.type = rec.info.get("SVTYPE")
        else:
            self.type = ""
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
        ide,
        cns,
        merged_vcf,
        ref_gene,
        out_dir,
    ):
        """Creation of sample object."""
        # pylint: disable=line-too-long
        self.id = ide
        self.cns = cns
        self.gen_reg = ref_gene
        self.out_dir = os.path.join(out_dir)
        self.vcf = str(merged_vcf)
        self.segs = os.path.join(self.out_dir, "%s_segs.csv" % self.id)
        self.circos_out = os.path.join(self.out_dir, "%s_circos.png" % self.id)
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
            self.out_dir, "%s_circos_svs.tsv" % self.id
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
        script_path = os.path.join(
             os.path.dirname(os.path.abspath(__file__)), "make_circos.r"
         )
        #script_path="/home/danielavt/toil_circosigv/toil_circosigv/make_circos.r"
        cmd = [
            script_path,
            circos_sv_file,
            self.id,
            self.gen_reg,
            self.segs,
            self.circos_out,
        ]
        print(" ".join(cmd)) 


if __name__ == "__main__":
    args = parser.parse_args()
    sampl= Sample(args.ide,args.cns,args.vcf,args.ref,args.out)
    sampl.runcircos()